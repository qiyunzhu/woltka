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
from tempfile import mkdtemp

from woltka.classify import (
    classify, assign, assign_none, assign_free, assign_rank, count, majority,
    write_profile, prep_table)


class ClassifyTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_classify(self):
        # simplest gotu workflow
        # for other tests see test_cli.classify
        input_path = join(self.datdir, 'align', 'bowtie2')
        output_path = join(self.tmpdir, 'tmp.tsv')
        classify(input_path, output_path)
        with open(output_path, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'bowtie2.gotu.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)
        remove(output_path)

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

    def test_write_profile(self):
        # default mode
        data = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                'S3': {'G2': 3, 'G5': 5}}
        fp = join(self.tmpdir, 'profile.tsv')
        with open(fp, 'w') as f:
            write_profile(f, data)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3',
               'G1\t4\t2\t0',
               'G2\t5\t0\t3',
               'G3\t8\t0\t0',
               'G4\t0\t3\t0',
               'G5\t0\t7\t5']
        self.assertListEqual(obs, exp)

        # with taxon names
        named = {'G1': 'Actinobacteria',
                 'G2': 'Firmicutes',
                 'G3': 'Bacteroidetes',
                 'G4': 'Cyanobacteria'}
        with open(fp, 'w') as f:
            write_profile(f, data, named=named)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3',
               'Actinobacteria\t4\t2\t0',
               'Firmicutes\t5\t0\t3',
               'Bacteroidetes\t8\t0\t0',
               'Cyanobacteria\t0\t3\t0',
               'G5\t0\t7\t5']
        self.assertListEqual(obs, exp)

        # with sample Ids
        samples = ['S3', 'S1']
        with open(fp, 'w') as f:
            write_profile(f, data, samples=samples)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS3\tS1',
               'G1\t0\t4',
               'G2\t3\t5',
               'G3\t0\t8',
               'G4\t0\t0',
               'G5\t5\t0']
        self.assertListEqual(obs, exp)
        remove(fp)

    def test_prep_table(self):
        # default mode
        data = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                'S3': {'G2': 3, 'G5': 5}}
        obs0, obs1, obs2 = prep_table(data)
        exp0 = [[4, 2, 0],
                [5, 0, 3],
                [8, 0, 0],
                [0, 3, 0],
                [0, 7, 5]]
        self.assertListEqual(obs0, exp0)
        self.assertListEqual(obs1, ['G1', 'G2', 'G3', 'G4', 'G5'])
        self.assertListEqual(obs2, ['S1', 'S2', 'S3'])

        # with sample Ids
        samples = ['S3', 'S1']
        obs0, _, obs2 = prep_table(data, samples=samples)
        exp0 = [[0, 4],
                [3, 5],
                [0, 8],
                [0, 0],
                [5, 0]]
        self.assertListEqual(obs0, exp0)
        self.assertListEqual(obs2, ['S3', 'S1'])


if __name__ == '__main__':
    main()
