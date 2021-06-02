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
    assign_none, assign_free, assign_rank, count_taxa, count_taxa_strat,
    total_taxa, majority)


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
        exp = ['G1', 'G2']
        self.assertListEqual(sorted(obs), exp)

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

        # assign one sub to nothing
        obs = assign_free({'Gx'}, **kwargs)
        self.assertIsNone(obs)

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
        exp = ['T1', 'T1', 'T2']
        self.assertListEqual(sorted(obs), exp)

        # assign above rank
        obs = assign_rank({'G1', 'G2', 'G3'}, 'general', tree=tree,
                          rankdic=rankdic, root=None, above=True)
        self.assertEqual(obs, 'T0')

        # assign above rank but there is a None
        obs = assign_rank({'G1', 'G2', 'G3', None}, 'general', tree=tree,
                          rankdic=rankdic, root=None, above=True)
        self.assertEqual(obs, None)

    def test_count_taxa(self):
        # G1, G4, G6: Ecoli
        # G2, G5: Cdiff
        # G3: Strep
        # G7: none
        subque = [['G1'],
                  ['G2', 'G3'],
                  ['G3', 'G4', 'G5', 'G6', 'G7'],
                  ['G4', 'G6'],
                  ['G7']]
        taxque = ['Ecoli',
                  ['Cdiff', 'Strep'],
                  ['Strep', 'Ecoli', 'Cdiff', 'Ecoli', None],
                  'Ecoli',
                  None]
        obs = count_taxa(subque, taxque)
        exp = {'Ecoli': {'G1': 1,   'G4': 0.75, 'G6': 0.75},
               'Cdiff': {'G2': 0.5, 'G5': 0.25},
               'Strep': {'G3': 0.75}}
        self.assertDictEqual(obs, exp)

    def test_count_taxa_strat(self):
        strata = {'seq1': 'Ecoli',
                  'seq2': 'Ecoli',
                  'seq3': 'Cdiff',
                  'seq4': 'Strep',
                  'seq5': 'Strep',
                  'seq6': 'Ecoli'}
        qryque = ['seq1',
                  'seq2',
                  'seq3',
                  'seq4',
                  'seq5',
                  'seq6',
                  'seq9']
        subque = [['G1'],
                  ['G2', 'G3'],
                  ['G3', 'G4', 'G5'],
                  ['G4', 'G6'],
                  ['G7'],
                  ['G3'],
                  ['G2', 'G6']]
        taxque = ['ligase',
                  ['polymerase', 'nuclease'],
                  ['nuclease', 'ligase', None],
                  'ligase',
                  None,
                  'nuclease',
                  ['polymerase', 'ligase']]
        obs = count_taxa_strat(qryque, subque, taxque, strata)
        exp = {('Ecoli', 'ligase'):     {'G1': 1},
               ('Ecoli', 'polymerase'): {'G2': 0.5},
               ('Ecoli', 'nuclease'):   {'G3': 1.5},
               ('Cdiff', 'nuclease'):   {'G3': 0.5},
               ('Cdiff', 'ligase'):     {'G4': 0.5},
               ('Strep', 'ligase'):     {'G4': 0.5, 'G6': 0.5}}
        self.assertDictEqual(obs, exp)

    def test_total_taxa(self):
        # without sizes
        counts = {'Ecoli': {'G1': 12, 'G2':  9, 'G3': 5},
                  'Cdiff': {'G4': 8,  'G5': 15}}
        obs = counts.copy()
        total_taxa(obs)
        exp = {'Ecoli': 26, 'Cdiff': 23}
        self.assertDictEqual(obs, exp)

        # with sizes
        sizes = {'G1': 100, 'G2': 80, 'G3': 75, 'G4': 120, 'G5': 50}
        obs = counts.copy()
        total_taxa(obs, sizes)
        exp = {'Ecoli': 0.29916667, 'Cdiff': 0.36666667}
        for taxon in exp:
            self.assertAlmostEqual(obs[taxon], exp[taxon])

        # subject not found
        del(sizes['G5'])
        obs = counts.copy()
        with self.assertRaises(ValueError) as ctx:
            total_taxa(counts, sizes)
        self.assertEqual(str(ctx.exception), (
            'One or more subjects are not found in the size map.'))

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
