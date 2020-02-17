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
from shutil import rmtree
from tempfile import mkdtemp

from woltka.util import readzip
from woltka.ordinal import match_read_gene, read_gene_coords, whether_prefix


class OrdinalTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

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
        obs = match_read_gene(queue, lens, th=0.8)
        exp = {'r1': {'g1'},
               'r5': {'g2'},
               'r6': {'g2'},
               'r8': {'g3'}}
        self.assertDictEqual(obs, exp)

        # threashold = 50%
        obs = match_read_gene(queue, lens, th=0.5)
        exp = {'r1': {'g1'},
               'r2': {'g1'},
               'r3': {'g1'},
               'r5': {'g2'},
               'r6': {'g2'},
               'r7': {'g3'},
               'r8': {'g3'},
               'r9': {'g3'}}
        self.assertDictEqual(obs, exp)

    def test_read_gene_coords(self):
        # simple case
        tbl = ('## GCF_000123456',
               '# NC_123456',
               '1	5	384',
               '2	410	933',
               '# NC_789012',
               '1	912	638',
               '2	529	75')
        obs = read_gene_coords(tbl, sort=True)
        exp = {'NC_123456': [
            (5,   True, True, '1'), (384, False, True, '1'),
            (410, True, True, '2'), (933, False, True, '2')],
               'NC_789012': [
            (75,  True, True, '2'), (529, False, True, '2'),
            (638, True, True, '1'), (912, False, True, '1')]}
        self.assertDictEqual(obs, exp)

        # real coords file
        fp = join(self.datdir, 'function', 'coords.txt.xz')
        with readzip(fp) as f:
            obs = read_gene_coords(f, sort=True)
        self.assertEqual(len(obs), 107)
        obs_ = obs['G000006745']
        self.assertEqual(len(obs_), 7188)
        self.assertTupleEqual(obs_[0], (372,  True,  True, '1'))
        self.assertTupleEqual(obs_[1], (806,  False, True, '1'))
        self.assertTupleEqual(obs_[2], (816,  True,  True, '2'))
        self.assertTupleEqual(obs_[3], (2177, False, True, '2'))

    def test_whether_prefix(self):
        # gene Ids are indices
        coords = {'NC_123456': [
            (5,   True, True, '1'), (384, False, True, '1'),
            (410, True, True, '2'), (933, False, True, '2')],
                  'NC_789012': [
            (75,  True, True, '2'), (529, False, True, '2'),
            (638, True, True, '1'), (912, False, True, '1')]}
        self.assertTrue(whether_prefix(coords))

        # gene Ids are unique accessions
        coords = {'NC_123456': [
            (5,   True, True,  'NP_135792.1'),
            (384, False, True, 'NP_135792.1'),
            (410, True, True,  'NP_246801.2'),
            (933, False, True, 'NP_246801.2')],
                  'NC_789012': [
            (75,  True, True,  'NP_258147.1'),
            (529, False, True, 'NP_258147.1'),
            (638, True, True,  'NP_369258.2'),
            (912, False, True, 'NP_369258.2')]}
        self.assertFalse(whether_prefix(coords))


if __name__ == '__main__':
    main()
