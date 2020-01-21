#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp, mkstemp

from woltk.ordinal import (
    match_read_gene, read_gene_table, readmap_to_profile)


class OrdinalTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_match_read_gene(self):
        # illustration of map

        #  reads:             ---------r1---------
        #                           ---------r2---------
        #                               ---------r3---------
        #                                 ---------r4---------
        #                                         ---------r5---------
        #                                         ---------r6---------
        # genome:  1 ====>>>>>>>>>>>g1>>>>>>>>>>>>===>>>>>>>>>>>>>>g2>> 50

        #  reads:             ---------r7---------
        #                          ---------r8---------
        #                                           ------r9------
        # genome: 51 >>>>>>>>>>>===>>>>>>>>>>>>>>g3>>>>>>>>>>>>>======= 100
        # --------------------------------------------------

        # gene table
        genes = [('g1',  5, 29),
                 ('g2', 33, 61),
                 ('g3', 65, 94)]
        # read map
        reads = [('r1', 10, 29),
                 ('r2', 16, 35),
                 ('r3', 20, 39),
                 ('r4', 22, 41),
                 ('r5', 30, 49),
                 ('r6', 30, 49),  # identical
                 ('r7', 60, 79),
                 ('r8', 65, 84),
                 ('r9', 82, 95)]  # shorter

        # length map
        lens = {'r{}'.format(i): 20 for i in range(1, 10)}

        # flatten lists
        genes = [x for id_, start, end in genes for x in
                 ((start, True, True, id_),
                  (end,  False, True, id_))]
        reads = [x for id_, start, end in reads for x in
                 ((start, True, False, id_),
                  (end,  False, False, id_))]

        queue = sorted(genes + reads, key=lambda x: x[0])

        # default (threshold = 80%)
        obs = match_read_gene(queue, lens)
        exp = {'g1': 1, 'g2': 2, 'g3': 1}
        self.assertDictEqual(obs, exp)

        # return read map instead of counts
        obs = match_read_gene(queue, lens, th=0.8, ismap=True)
        exp = {'r1': {'g1'},
               'r5': {'g2'},
               'r6': {'g2'},
               'r8': {'g3'}}
        self.assertDictEqual(obs, exp)

        # threashold = 50%
        obs = match_read_gene(queue, lens, th=0.5, ismap=True)
        exp = {'r1': {'g1'},
               'r2': {'g1'},
               'r3': {'g1'},
               'r5': {'g2'},
               'r6': {'g2'},
               'r7': {'g3'},
               'r8': {'g3'},
               'r9': {'g3'}}
        self.assertDictEqual(obs, exp)

    def test_read_gene_table(self):
        tbl = ('## GCF_000123456',
               '# NC_123456',
               '1	5	384',
               '2	410	933',
               '# NC_789012',
               '1	912	638',
               '2	529	75')
        obs = read_gene_table(tbl)
        exp = {'NC_123456': [
            (5,   True, True, '1'), (384, False, True, '1'),
            (410, True, True, '2'), (933, False, True, '2')],
               'NC_789012': [
            (638, True, True, '1'), (912, False, True, '1'),
            (75,  True, True, '2'), (529, False, True, '2')]}
        self.assertDictEqual(obs, exp)

    def test_readmap_to_profile(self):
        # no ambiguity
        rids = ['R1', 'R2', 'R3', 'R4']
        profile = {'G1': [0, 1, 2], 'G2': [1, 3]}
        readmap_to_profile(profile, rids, True)
        self.assertDictEqual(profile, {'G1': 3, 'G2': 2})

        # R1 occurs 3 times; R2 occurs 2 times
        rids = ['R1', 'R2', 'R3', 'R1', 'R1', 'R2']

        # drop R1 and R2
        profile = {'G1': [0, 1], 'G2': [1, 2, 3]}
        readmap_to_profile(profile, rids, False)
        self.assertDictEqual(profile, {'G2': 1})

        # normalize R1 by 3, R2 by 2
        profile = {'G1': [0, 1], 'G2': [1, 2, 3]}
        readmap_to_profile(profile, rids, True)
        self.assertDictEqual(profile, {'G1': 1, 'G2': 2})


if __name__ == '__main__':
    main()
