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
    assign, count, majority,  write_profile, prep_table)


class ClassifyTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_classify(self):
        # TODO
        pass

    def test_assign(self):
        # without a tree
        # assign to self
        self.assertEqual(assign({'G1'}), 'G1')

        # count subjects
        obs = assign({'G1', 'G2'}, ambig=True)
        exp = {'G1': 1, 'G2': 1}
        self.assertDictEqual(obs, exp)

        # cannot assign
        obs = assign({'G1', 'G2', 'G3'})
        self.assertIsNone(obs)

        # rank has no effect if no tree
        obs = assign({'G1', 'G2', 'G3'}, rank='hello')
        self.assertIsNone(obs)

        # with a tree
        subs = {'G1', 'G2', 'G3'}
        tree = {'G1': 'T1', 'G2': 'T1', 'G3': 'T2',
                'T1': 'T0', 'T2': 'T0', 'T0': 'T0'}
        rankd = {'T1': 'general', 'T2': 'general', 'T0': 'marshal'}

        # free-rank assign
        obs = assign(subs, 'free', tree=tree)
        self.assertEqual(obs, 'T0')

        # tree not used
        obs = assign(subs, 'none', tree=tree)
        self.assertIsNone(obs)

        # fixed rank assign
        kwargs = {'tree': tree, 'rankd': rankd}
        obs = assign(subs, 'marshal', **kwargs)
        self.assertEqual(obs, 'T0')
        obs = assign(subs, 'general', **kwargs)
        self.assertIsNone(obs)
        obs = assign({'G1', 'G2'}, 'general', **kwargs)
        self.assertEqual(obs, 'T1')

        # rank that does not exist
        obs = assign(subs, 'admiral', **kwargs)
        self.assertIsNone(obs)

        # TODO: add more tests

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
