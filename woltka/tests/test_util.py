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

from woltka.util import (
    add_dict, update_dict, sum_dict, intize, intize_list, intize_dict,
    delnone, allkeys, count_list, last_value, feature_count)


class UtilTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_add_dict(self):
        d = {'a': 1, 'b': 2}
        add_dict(d, 'c', 3)
        self.assertDictEqual(d, {'a': 1, 'b': 2, 'c': 3})
        add_dict(d, 'b', 2)
        self.assertDictEqual(d, {'a': 1, 'b': 2, 'c': 3})
        with self.assertRaises(AssertionError) as ctx:
            add_dict(d, 'a', 0)
        self.assertEqual(str(ctx.exception), (
            'Conflicting values found for "a".'))

    def test_update_dict(self):
        d0 = {'a': 1, 'b': 2}
        d1 = {'c': 3, 'd': 4}
        update_dict(d0, d1)
        self.assertDictEqual(d0, {'a': 1, 'b': 2, 'c': 3, 'd': 4})
        d2 = {'e': 5, 'a': 1}
        update_dict(d0, d2)
        self.assertDictEqual(d0, {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5})
        d3 = {'g': 6, 'b': 100}
        with self.assertRaises(AssertionError) as ctx:
            update_dict(d0, d3)
        self.assertEqual(str(ctx.exception), (
            'Conflicting values found for "b".'))

    def test_sum_dict(self):
        d0 = {'a': 1, 'b': 2}
        d1 = {'c': 3, 'd': 4}
        sum_dict(d0, d1)
        self.assertDictEqual(d0, {'a': 1, 'b': 2, 'c': 3, 'd': 4})
        d2 = {'a': 1}
        sum_dict(d0, d2)
        self.assertDictEqual(d0, {'a': 2, 'b': 2, 'c': 3, 'd': 4})
        d3 = {'e': 5, 'b': 3}
        sum_dict(d0, d3)
        self.assertDictEqual(d0, {'a': 2, 'b': 5, 'c': 3, 'd': 4, 'e': 5})

    def test_intize(self):
        self.assertEqual(intize(1), 1)
        self.assertEqual(intize(1.0), 1)
        self.assertEqual(intize(1.1), 1)
        self.assertEqual(intize(1.49), 1)
        self.assertEqual(intize(1.49999), 1)
        self.assertEqual(intize(1.49999999999), 2)
        self.assertEqual(intize(1.5), 2)
        self.assertEqual(intize(1.8), 2)
        self.assertEqual(intize(2.2), 2)
        self.assertEqual(intize(2.5), 2)
        self.assertEqual(intize(2.50001), 3)
        self.assertEqual(intize(2.50000000001), 2)
        self.assertEqual(intize(2.7), 3)
        self.assertEqual(intize(-1.0), -1)
        self.assertEqual(intize(-1.5), -2)
        self.assertEqual(intize(-1.49999999999), -2)

    def test_intize_list(self):
        lst = [0.1, 1.2, 2.4, 3.8, 4.5, 5.5]
        intize_list(lst)
        exp = [0, 1, 2, 4, 4, 6]
        self.assertListEqual(lst, exp)

    def test_intize_dict(self):
        dic = {'a': 1.0, 'b': 2.2, 'c': 3.6}
        intize_dict(dic)
        exp = {'a': 1, 'b': 2, 'c': 4}
        self.assertDictEqual(dic, exp)
        dic = {'a': -0.2, 'b': -3.3, 'c': 1.8, 'd': 0.4}
        intize_dict(dic)
        exp = {'b': -3, 'c': 2}
        self.assertDictEqual(dic, exp)
        dic = {'a': -0.2, 'b': -3.3, 'c': 1.8, 'd': 0.4}
        intize_dict(dic, zero=True)
        exp = {'a': 0, 'b': -3, 'c': 2, 'd': 0}
        self.assertDictEqual(dic, exp)
        dic = {'a': 2.5, 'b': 3.5, 'c': 1.49999999999, 'd': 1.50000000001}
        intize_dict(dic)
        exp = {'a': 2, 'b': 4, 'c': 2, 'd': 2}
        self.assertDictEqual(dic, exp)

    def test_delnone(self):
        dic = {'a': 1, 'b': 2, 'c': 3}
        obs = dic.copy()
        delnone(obs)
        self.assertDictEqual(obs, dic)
        obs[None] = 4
        delnone(obs)
        self.assertDictEqual(obs, dic)

    def test_allkeys(self):
        data = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                'S3': {'G2': 3, 'G5': 5}}
        obs = allkeys(data)
        exp = {'G1', 'G2', 'G3', 'G4', 'G5'}
        self.assertSetEqual(obs, exp)

    def test_count_list(self):
        lst = [1, 1, 2, 3, 1, 2, 4, 3, 1]
        obs = count_list(lst)
        exp = {1: 4, 2: 2, 3: 2, 4: 1}
        self.assertDictEqual(obs, exp)

    def test_last_value(self):
        lst = ['a', 1, None, 'b', 2, None]
        self.assertEqual(last_value(lst), 2)
        self.assertIsNone(last_value([None, None]))

    def test_feature_count(self):
        self.assertTupleEqual(feature_count('A'),   ('A', 1))
        self.assertTupleEqual(feature_count('A:1'), ('A', 1))
        self.assertTupleEqual(feature_count('A:2'), ('A', 2))
        self.assertTupleEqual(feature_count('A:0'), ('A:0', 1))
        self.assertTupleEqual(feature_count('A:x'), ('A:x', 1))
        self.assertTupleEqual(feature_count(':1'),  (':1', 1))
        self.assertTupleEqual(feature_count('::'),  ('::', 1))
        self.assertTupleEqual(feature_count(''),    ('', 1))


if __name__ == '__main__':
    main()
