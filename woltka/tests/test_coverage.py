#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join, dirname, realpath
from os import makedirs
from shutil import rmtree
from tempfile import mkdtemp

from woltka.coverage import (
    merge_ranges, parse_ranges, calc_coverage, write_coverage)


class CoverageTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_merge_ranges(self):
        # ranges that are sorted and overlapped
        obs = merge_ranges([1, 3, 2, 4, 6, 8, 7, 9])
        exp = [1, 4, 6, 9]
        self.assertListEqual(obs, exp)

        # unsorted ranges
        obs = merge_ranges([4, 6, 1, 4, 5, 9])
        exp = [1, 9]
        self.assertListEqual(obs, exp)

        # ranges that are connected (but not overlapped)
        obs = merge_ranges([1, 2, 3, 4, 5, 6])
        exp = [1, 6]
        self.assertListEqual(obs, exp)

        # empty range list
        self.assertListEqual(merge_ranges([]), [])

    def test_parse_ranges(self):
        rmap = {'S1': (['R1', 'R2', 'R3'], [
                    {'G1': [1, 100]},
                    {'G1': [1, 50, 251, 300]},
                    {'G1': [1, 50], 'G2': [101, 150]}]),
                'S2': (['R1', 'R3', 'R4'], [
                    {'G2': [151, 200, 51, 100]},
                    {'G3': [76, 125], 'G4': [26, 75, 101, 150]},
                    {'G2': [26, 75], 'G3': [1, 50]}]),
                'S3': (['R2', 'R4', 'R5'], [
                    {'G1': [1, 50], 'G2': [51, 100]},
                    {'G1': [51, 100, 151, 200]},
                    {'G1': [101, 150], 'G2': [1, 50]}])}

        # add ranges but don't merge
        obs = {}
        parse_ranges(rmap, obs)
        exp = {('S1', 'G1'): [8, [1, 100, 1, 50, 251, 300, 1, 50]],
               ('S1', 'G2'): [2, [101, 150]],
               ('S2', 'G2'): [6, [151, 200, 51, 100, 26, 75]],
               ('S2', 'G3'): [4, [76, 125, 1, 50]],
               ('S2', 'G4'): [4, [26, 75, 101, 150]],
               ('S3', 'G1'): [8, [1, 50, 51, 100, 151, 200, 101, 150]],
               ('S3', 'G2'): [4, [51, 100, 1, 50]]}
        self.assertDictEqual(obs, exp)

        # merge while adding
        obs = {}
        parse_ranges(rmap, obs, chunk=1)
        exp = {('S1', 'G1'): [0, [1, 100, 251, 300]],
               ('S1', 'G2'): [0, [101, 150]],
               ('S2', 'G2'): [0, [26, 100, 151, 200]],
               ('S2', 'G3'): [0, [1, 50, 76, 125]],
               ('S2', 'G4'): [0, [26, 75, 101, 150]],
               ('S3', 'G1'): [0, [1, 200]],
               ('S3', 'G2'): [0, [1, 100]]}
        self.assertDictEqual(obs, exp)

    def test_calc_coverage(self):
        cov = {('S1', 'G1'): [8, [1, 100, 1, 50, 251, 300, 1, 50]],
               ('S1', 'G2'): [2, [101, 150]],
               ('S2', 'G2'): [6, [151, 200, 51, 100, 26, 75]],
               ('S2', 'G3'): [4, [76, 125, 1, 50]],
               ('S2', 'G4'): [4, [26, 75, 101, 150]],
               ('S3', 'G1'): [8, [1, 50, 51, 100, 151, 200, 101, 150]],
               ('S3', 'G2'): [4, [51, 100, 1, 50]]}
        obs = calc_coverage(cov)
        exp = {'S1': {'G1': [1, 100, 251, 300],
                      'G2': [101, 150]},
               'S2': {'G2': [26, 100, 151, 200],
                      'G3': [1, 50, 76, 125],
                      'G4': [26, 75, 101, 150]},
               'S3': {'G1': [1, 200],
                      'G2': [1, 100]}}
        self.assertDictEqual(obs, exp)

    def test_write_coverage(self):
        cov = {'S1': {'G1': [1, 100, 251, 300],
                      'G2': [101, 150]},
               'S2': {'G2': [26, 100, 151, 200],
                      'G3': [1, 50, 76, 125],
                      'G4': [26, 75, 101, 150]},
               'S3': {'G1': [1, 200],
                      'G2': [1, 100]}}
        outdir = join(self.tmpdir, 'outdir')
        makedirs(outdir)
        write_coverage(cov, outdir)
        with open(join(outdir, 'S1.cov'), 'r') as f:
            obs = f.read().splitlines()
        exp = ['G1\t1\t100', 'G1\t251\t300', 'G2\t101\t150']
        self.assertListEqual(obs, exp)
        with open(join(outdir, 'S2.cov'), 'r') as f:
            obs = f.read().splitlines()
        exp = ['G2\t26\t100', 'G2\t151\t200', 'G3\t1\t50',
               'G3\t76\t125', 'G4\t26\t75', 'G4\t101\t150']
        self.assertListEqual(obs, exp)
        with open(join(outdir, 'S3.cov'), 'r') as f:
            obs = f.read().splitlines()
        exp = ['G1\t1\t200', 'G2\t1\t100']
        self.assertListEqual(obs, exp)
        rmtree(outdir)


if __name__ == '__main__':
    main()
