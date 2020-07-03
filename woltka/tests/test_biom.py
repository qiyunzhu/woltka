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
from io import StringIO

import numpy as np
import pandas as pd
from biom import load_table, Table
from pandas.testing import assert_frame_equal

from woltka.biom import table_to_biom, biom_to_table, write_biom, filter_biom
from woltka.table import prep_table


class BiomTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_table_to_biom(self):
        data = [[4, 2, 0],
                [5, 0, 3],
                [8, 0, 0],
                [0, 3, 0],
                [0, 7, 5]]
        observs = ['G1', 'G2', 'G3', 'G4', 'G5']
        samples = ['S1', 'S2', 'S3']
        metadata = [
            {'Name': 'Actinobacteria', 'Rank': 'phylum', 'Lineage': '2;72;74'},
            {'Name': 'Firmicutes',     'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': 'Bacteroidetes',  'Rank': 'phylum', 'Lineage': '2;70'},
            {'Name': 'Cyanobacteria',  'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': '',               'Rank': '',       'Lineage': ''}]
        obs = table_to_biom(data, observs, samples, metadata)
        exp = pd.read_csv(StringIO(
            '\tS1\tS2\tS3\tName\tRank\tLineage\n'
            'G1\t4\t2\t0\tActinobacteria\tphylum\t2;72;74\n'
            'G2\t5\t0\t3\tFirmicutes\tphylum\t2;72\n'
            'G3\t8\t0\t0\tBacteroidetes\tphylum\t2;70\n'
            'G4\t0\t3\t0\tCyanobacteria\tphylum\t2;72\n'
            'G5\t0\t7\t5\t\t\t\n'),
            sep='\t', index_col=0, na_filter=False)
        assert_frame_equal(
            obs.to_dataframe(dense=True).astype(int), exp.iloc[:, :3])
        assert_frame_equal(
            obs.metadata_to_dataframe('observation')[[
                'Name', 'Rank', 'Lineage']], exp.iloc[:, -3:])

    def test_biom_to_table(self):
        data = [[4, 2, 0],
                [5, 0, 3],
                [8, 0, 0],
                [0, 3, 0],
                [0, 7, 5]]
        observs = ['G1', 'G2', 'G3', 'G4', 'G5']
        samples = ['S1', 'S2', 'S3']
        metadata = [
            {'Name': 'Actinobacteria', 'Rank': 'phylum', 'Lineage': '2;72;74'},
            {'Name': 'Firmicutes',     'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': 'Bacteroidetes',  'Rank': 'phylum', 'Lineage': '2;70'},
            {'Name': 'Cyanobacteria',  'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': '',               'Rank': '',       'Lineage': ''}]
        table = table_to_biom(data, observs, samples, metadata)
        obs = biom_to_table(table)
        self.assertListEqual(obs[0], data)
        self.assertListEqual(obs[1], observs)
        self.assertListEqual(obs[2], samples)
        self.assertListEqual(obs[3], metadata)

    def test_write_biom(self):
        profile = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                   'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                   'S3': {'G2': 3, 'G5': 5}}
        exp = table_to_biom(*prep_table(profile))
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
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        obs = filter_biom(table, th=0.5)
        exp = Table(*map(np.array, prep_table({
            'S1': {},
            'S2': {'G5': 7},
            'S3': {'G5': 5}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # empty BIOM table cannot be directly compared
        obs = filter_biom(table, th=10)
        self.assertTupleEqual(obs.to_dataframe(True).shape, (0, 3))


if __name__ == '__main__':
    main()
