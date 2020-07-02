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

import numpy as np
import pandas as pd
from biom import load_table, Table
from numpy.testing import assert_array_equal
from pandas.testing import assert_frame_equal

from woltka.biom import profile_to_biom, write_biom, filter_biom
from woltka.file import prep_table


class BiomTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_profile_to_biom(self):
        # default mode
        profile = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                   'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                   'S3': {'G2': 3, 'G5': 5}}
        table = profile_to_biom(profile)
        samples = ['S1', 'S2', 'S3']
        assert_array_equal(table.ids('sample'), samples)
        features = ['G1', 'G2', 'G3', 'G4', 'G5']
        assert_array_equal(table.ids('observation'), features)
        data = [[4, 2, 0], [5, 0, 3], [8, 0, 0], [0, 3, 0], [0, 7, 5]]
        obs = table.to_dataframe(dense=True).astype(int).values
        assert_array_equal(obs, data)

        # with sample Ids
        table = profile_to_biom(profile, samples=['S3', 'S1'])
        obs = table.to_dataframe(dense=True).astype(int)
        exp = pd.DataFrame([[0, 4], [3, 5], [0, 8], [0, 0], [5, 0]],
                           index=features, columns=['S3', 'S1'])
        assert_frame_equal(obs, exp)

        # some sample Ids are not in profile
        table = profile_to_biom(profile, samples=['S3', 'S0', 'S1'])
        obs = table.to_dataframe(dense=True).astype(int)
        assert_frame_equal(obs, exp)

        # with taxon names
        namedic = {'G1': 'Actinobacteria',
                   'G2': 'Firmicutes',
                   'G3': 'Bacteroidetes',
                   'G4': 'Cyanobacteria'}
        table = profile_to_biom(profile, namedic=namedic)
        obs = table.to_dataframe(dense=True).astype(int)
        exp = pd.DataFrame(data, features, samples)
        assert_frame_equal(obs, exp)
        obs = table.metadata_to_dataframe('observation')['Name']
        names = ['Actinobacteria', 'Firmicutes', 'Bacteroidetes',
                 'Cyanobacteria', None]
        assert_array_equal(obs, names)

        # with taxon names to replace Ids
        table = profile_to_biom(profile, namedic=namedic, name_as_id=True)
        obs = table.to_dataframe(dense=True).astype(int)
        exp = pd.DataFrame(data, names[:4] + ['G5'], samples)
        assert_frame_equal(obs, exp)

        # with ranks
        rankdic = {'G1': 'class', 'G2': 'phylum', 'G4': 'phylum'}
        table = profile_to_biom(profile, rankdic=rankdic)
        obs = table.metadata_to_dataframe('observation')['Rank']
        exp = ['class', 'phylum', None, 'phylum', None]
        assert_array_equal(obs, exp)

        # with lineages
        tree = {'G1': '74', '74': '72', 'G2': '72', 'G3': '70', 'G4': '72',
                'G5':  '1', '72':  '2', '70':  '2',  '2':  '1',  '1':  '1'}
        table = profile_to_biom(profile, tree=tree)
        obs = table.metadata_to_dataframe('observation')['Lineage']
        exp = ['2;72;74', '2;72', '2;70', '2;72', None]
        assert_array_equal(obs, exp)

        # with lineages and names as Ids
        namedic.update({
            '74': 'Actino', '72': 'Terra', '70': 'FCB', '2': 'Bacteria'})
        table = profile_to_biom(
            profile, tree=tree, namedic=namedic, name_as_id=True)
        obs = table.metadata_to_dataframe('observation')['Lineage']
        exp = ['Bacteria;Terra;Actino', 'Bacteria;Terra', 'Bacteria;FCB',
               'Bacteria;Terra', None]
        assert_array_equal(obs, exp)

        # with stratification
        profile = {'S1': {('A', 'G1'): 4,
                          ('A', 'G2'): 5,
                          ('B', 'G1'): 8},
                   'S2': {('A', 'G1'): 2,
                          ('B', 'G1'): 3,
                          ('B', 'G2'): 7},
                   'S3': {('B', 'G3'): 3,
                          ('C', 'G2'): 5}}
        table = profile_to_biom(profile)
        obs = table.to_dataframe(dense=True).astype(int)
        data = [[4, 2, 0], [5, 0, 0], [8, 3, 0], [0, 7, 0], [0, 0, 3],
                [0, 0, 5]]
        index = ['A|G1', 'A|G2', 'B|G1', 'B|G2', 'B|G3', 'C|G2']
        exp = pd.DataFrame(data, index=index, columns=samples)
        assert_frame_equal(obs, exp)

    def test_write_biom(self):
        profile = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                   'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                   'S3': {'G2': 3, 'G5': 5}}
        exp = profile_to_biom(profile)
        fp = join(self.tmpdir, 'tmp.biom')
        write_biom(exp, fp)
        obs = load_table(fp)
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')
        remove(fp)

    def test_filter_biom(self):
        table = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8},
            'S2': {'G1': 2, 'G4': 3, 'G5': 7},
            'S3': {'G2': 3, 'G5': 5}})))
        obs = filter_biom(table, th=3)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8},
            'S2': {'G4': 3, 'G5': 7},
            'S3': {'G2': 3, 'G5': 5}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        obs = filter_biom(table, th=4)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8},
            'S2': {'G5': 7},
            'S3': {'G5': 5}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        obs = filter_biom(table, th=6)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G3': 8},
            'S2': {'G5': 7},
            'S3': {}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        obs = filter_biom(table, th=0.25)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G2': 5, 'G3': 8},
            'S2': {'G4': 3, 'G5': 7},
            'S3': {'G2': 3, 'G5': 5}})))
        print(obs.to_dataframe(dense=True))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        obs = filter_biom(table, th=0.5)
        exp = Table(*map(np.array, prep_table({
            'S1': {},
            'S2': {'G5': 7},
            'S3': {'G5': 5}})))
        print(obs.to_dataframe(dense=True))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # empty BIOM table cannot be directly compared
        obs = filter_biom(table, th=10)
        self.assertTupleEqual(obs.to_dataframe(True).shape, (0, 3))


if __name__ == '__main__':
    main()
