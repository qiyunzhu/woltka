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

from woltka.classify import (
    assign, assign_none, assign_free, assign_rank, count, majority,
    strip_index, demultiplex)


class ClassifyTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_assign(self):
        # count subjects
        obs = assign({'G1', 'G2', 'G3'}, rank='none', ambig=True)
        self.assertDictEqual(obs, {'G1': 1, 'G2': 1, 'G3': 1})

        # free-rank assign
        tree = {'G1': 'T1', 'G2': 'T1', 'G3': 'T2',
                'T1': 'T0', 'T2': 'T0', 'T0': 'T0'}
        obs = assign({'G1', 'G2', 'G3'}, rank='free', tree=tree)
        self.assertEqual(obs, 'T0')
        obs = assign({'G1', 'G2', 'G3'}, rank='free', tree=tree, root='T0')
        self.assertIsNone(obs)

        # fixed rank assign
        rankd = {'T1': 'general', 'T2': 'general', 'T0': 'marshal'}
        kwargs = {'tree': tree, 'rankd': rankd, 'root': 'T0'}
        obs = assign({'G1', 'G2', 'G3'}, 'marshal', **kwargs)
        self.assertEqual(obs, 'T0')
        obs = assign({'G1', 'G2', 'G3'}, 'general', **kwargs)
        self.assertIsNone(obs)
        obs = assign({'G1', 'G2'}, 'general', **kwargs)
        self.assertEqual(obs, 'T1')

    def test_assign_none(self):
        # assign to self
        self.assertEqual(assign_none({'G1'}), 'G1')

        # count subjects
        obs = assign_none({'G1', 'G2'}, ambig=True)
        exp = {'G1': 1, 'G2': 1}
        self.assertDictEqual(obs, exp)

        # cannot assign
        obs = assign_none({'G1', 'G2', 'G3'})
        self.assertIsNone(obs)

    def test_assign_free(self):
        tree = {'G1': 'T1', 'G2': 'T1', 'G3': 'T2',
                'T1': 'T0', 'T2': 'T0', 'T0': 'T0'}
        kwargs = {'tree': tree, 'root': 'T0'}

        # free-rank assignment
        obs = assign_free({'G1', 'G2'}, **kwargs)
        self.assertEqual(obs, 'T1')

        # root -> None
        obs = assign_free({'G1', 'G2', 'G3'}, **kwargs)
        self.assertIsNone(obs)

        # assign one sub to a taxon
        obs = assign_free({'G1'}, **kwargs)
        self.assertEqual(obs, 'T1')

        # assign one sub to itself
        obs = assign_free({'G1'}, **kwargs, subok=True)
        self.assertEqual(obs, 'G1')

    def test_assign_rank(self):
        tree = {'G1': 'T1', 'G2': 'T1', 'G3': 'T2',
                'T1': 'T0', 'T2': 'T0', 'T0': 'T0'}
        rankd = {'T1': 'general', 'T2': 'general', 'T0': 'marshal'}
        kwargs = {'tree': tree, 'rankd': rankd, 'root': 'T0'}

        # fixed-rank assignment
        obs = assign_rank({'G1', 'G2', 'G3'}, 'marshal', **kwargs)
        self.assertEqual(obs, 'T0')
        obs = assign_rank({'G1', 'G2'}, 'general', **kwargs)
        self.assertEqual(obs, 'T1')
        obs = assign_rank({'G1'}, 'general', **kwargs)
        self.assertEqual(obs, 'T1')

        # rank that does not exist
        obs = assign_rank({'G1', 'G2', 'G3'}, 'admiral', **kwargs)
        self.assertIsNone(obs)

        # cannot find consensus at rank
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', **kwargs)
        self.assertIsNone(obs)

        # majority rule
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', **kwargs, major=0.6)
        self.assertEqual(obs, 'T1')

        # consider ambiguity
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', **kwargs, ambig=True)
        exp = {'T1': 2, 'T2': 1}
        self.assertDictEqual(obs, exp)

        # assign above rank
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', tree=tree,
                          rankd=rankd, root=None, above=True)
        self.assertEqual(obs, 'T0')

    def test_count(self):
        # unique match
        matches = {'seq1': 'a', 'seq2': 'a',  'seq3': 'b',
                   'seq4': 'c', 'seq5': None, 'seq6': 'b'}
        obs = count(matches)
        exp = {'a': 2, 'b': 2, 'c': 1, None: 1}
        self.assertDictEqual(obs, exp)

        # multiple matches
        matches = {'seq1': {'a': 1, 'b': 2, 'c': 3},
                   'seq2': {'a': 2, 'b': 5},
                   'seq3': {'d': 4}}
        obs = count(matches)
        exp = {'a': 1 / 6 + 2 / 7,
               'b': 2 / 6 + 5 / 7,
               'c': 3 / 6,
               'd': 4 / 4}
        for key in obs:
            self.assertAlmostEqual(obs[key], exp[key])

    def test_majority(self):
        obs = majority([1, 1, 1, 2, 2], th=0.6)
        self.assertEqual(obs, 1)
        obs = majority([1, 1, 1, 2, 2], th=0.8)
        self.assertIsNone(obs)
        obs = majority([1, 2, 3])
        self.assertIsNone(obs)

    def test_strip_index(self):
        dic = {'R1': ['G1_1', 'G1_2', 'G2_3', 'G3'],
               'R2': ['G1_1', 'G1.3', 'G4_5', 'G4_x']}
        strip_index(dic)
        self.assertDictEqual(dic, {
            'R1': ['G1', 'G1',   'G2', 'G3'],
            'R2': ['G1', 'G1.3', 'G4', 'G4']})

    def test_demultiplex(self):
        # simple case
        dic = {'S1_R1': 5,
               'S1_R2': 12,
               'S1_R3': 3,
               'S2_R1': 10,
               'S2_R2': 8,
               'S2_R4': 7,
               'S3_R2': 15,
               'S3_R3': 1,
               'S3_R4': 5}
        obs = demultiplex(dic)
        exp = {'S1': {'R1': 5, 'R2': 12, 'R3': 3},
               'S2': {'R1': 10, 'R2': 8, 'R4': 7},
               'S3': {'R2': 15, 'R3': 1, 'R4': 5}}
        self.assertDictEqual(obs, exp)

        # change separator, no result
        obs = demultiplex(dic, sep='.')
        self.assertDictEqual(obs, {'': dic})

        # enforce sample Ids
        obs = demultiplex(dic, samples=['S1', 'S2', 'SX'])
        exp = {x: exp[x] for x in ['S1', 'S2']}
        self.assertDictEqual(obs, exp)


if __name__ == '__main__':
    main()
