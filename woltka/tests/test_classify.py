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
    assign_none, assign_free, assign_rank, count, count_strata, majority)


class ClassifyTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_assign_none(self):
        # assign to self
        self.assertEqual(assign_none({'G1'}), 'G1')

        # count subjects
        obs = assign_none({'G1', 'G2'})
        exp = {'G1': 1, 'G2': 1}
        self.assertDictEqual(obs, exp)

        # cannot assign
        obs = assign_none({'G1', 'G2', 'G3'}, uniq=True)
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
        rankdic = {'T1': 'general', 'T2': 'general', 'T0': 'marshal'}
        kwargs = {'tree': tree, 'rankdic': rankdic, 'root': 'T0'}

        # fixed-rank assignment
        obs = assign_rank({'G1', 'G2', 'G3'}, 'marshal', **kwargs)
        self.assertEqual(obs, 'T0')
        obs = assign_rank({'G1', 'G2'}, 'general', **kwargs)
        self.assertEqual(obs, 'T1')
        obs = assign_rank({'G1'}, 'general', **kwargs)
        self.assertEqual(obs, 'T1')

        # rank that does not exist
        obs = assign_rank({'G1', 'G2', 'G3'}, 'admiral', **kwargs, uniq=True)
        self.assertIsNone(obs)

        # cannot find consensus at rank
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', **kwargs, uniq=True)
        self.assertIsNone(obs)

        # majority rule
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', **kwargs, major=0.6)
        self.assertEqual(obs, 'T1')

        # consider ambiguity
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', **kwargs)
        exp = {'T1': 2, 'T2': 1}
        self.assertDictEqual(obs, exp)

        # assign above rank
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', tree=tree,
                          rankdic=rankdic, root=None, above=True)
        self.assertEqual(obs, 'T0')

        # assign above rank but there is a None
        obs = assign_rank({'G1', 'G2', 'G3', None}, 'general', tree=tree,
                          rankdic=rankdic, root=None, above=True)
        self.assertEqual(obs, None)

    def test_count(self):
        # unique match
        taxque = ['a', 'a', 'b', 'c', None, 'b', None]
        obs = count(taxque)
        exp = {'a': 2, 'b': 2, 'c': 1, None: 2}
        self.assertDictEqual(obs, exp)

        # multiple matches
        taxque = [{'a': 1, 'b': 2, 'c': 3},
                  {'a': 2, 'b': 5},
                  {'d': 4}]
        obs = count(taxque)
        exp = {'a': 1 / 6 + 2 / 7,
               'b': 2 / 6 + 5 / 7,
               'c': 3 / 6,
               'd': 4 / 4}
        for key in obs:
            self.assertAlmostEqual(obs[key], exp[key])

    def test_count_strata(self):
        # unique match
        strata = {'seq1': 'Ecoli',
                  'seq2': 'Ecoli',
                  'seq3': 'Cdiff',
                  'seq4': 'Strep',
                  'seq5': 'Strep',
                  'seq6': 'Ecoli'}
        matches = {'seq1': 'ligase',
                   'seq2': 'polymerase',
                   'seq3': 'nuclease',
                   'seq4': 'ligase',
                   'seq5': 'nuclease',
                   'seq6': 'ligase',
                   'seq0': 'nothing'}
        obs = count_strata(matches.keys(), matches.values(), strata)
        exp = {('Ecoli',     'ligase'): 2,
               ('Ecoli', 'polymerase'): 1,
               ('Cdiff',   'nuclease'): 1,
               ('Strep',     'ligase'): 1,
               ('Strep',   'nuclease'): 1}
        self.assertDictEqual(obs, exp)

        # multiple matches
        strata['seq7'] = 'Ecoli'
        matches['seq7'] = {'polymerase': 3, 'ligase': 1}
        obs = count_strata(matches.keys(), matches.values(), strata)
        exp = {('Ecoli',     'ligase'): 2.25,
               ('Ecoli', 'polymerase'): 1.75,
               ('Cdiff',   'nuclease'): 1,
               ('Strep',     'ligase'): 1,
               ('Strep',   'nuclease'): 1}
        self.assertDictEqual(obs, exp)

    def test_majority(self):
        obs = majority([1, 1, 1, 2, 2], th=0.6)
        self.assertEqual(obs, 1)
        obs = majority([1, 1, 1, 2, 2], th=0.8)
        self.assertIsNone(obs)
        obs = majority([1, 2, 3])
        self.assertIsNone(obs)
        obs = majority([])
        self.assertIsNone(obs)


if __name__ == '__main__':
    main()
