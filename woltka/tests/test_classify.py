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

from woltka.classify import count, majority, assign, prep_table


class ClassifyTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

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

    def test_assign(self):
        obs = assign(['g1', 'g2', 'g3'])
        self.assertIsNone(obs)

    def test_prep_table(self):
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


if __name__ == '__main__':
    main()
