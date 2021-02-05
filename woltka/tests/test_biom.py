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

from woltka.biom import (
    table_to_biom, biom_to_table, write_biom, filter_biom, round_biom,
    biom_add_metacol, collapse_biom)
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

    def test_round_biom(self):
        obs = Table(*map(np.array, prep_table({
            'S1': {'G1': 0.5, 'G2': 0.0, 'G3': 2.3, 'G4': 0.50000000001},
            'S2': {'G1': 1.5, 'G2': 0.2, 'G3': 1.49999999999, 'G4': 0.2},
            'S3': {'G1': 2.5, 'G2': 0.3, 'G3': 3.8, 'G4': 0.1}})))
        round_biom(obs)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 0, 'G3': 2},
            'S2': {'G1': 2, 'G3': 2},
            'S3': {'G1': 2, 'G3': 4}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

    def test_biom_add_metacol(self):
        obs = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8, 'G4': 0, 'G5': 3},
            'S2': {'G1': 1, 'G2': 8, 'G3': 0, 'G4': 7, 'G5': 4},
            'S3': {'G1': 0, 'G2': 2, 'G3': 3, 'G4': 5, 'G5': 0}})))
        self.assertIsNone(obs.metadata(axis='observation'))
        rankdic = {'G1': 'S', 'G2': 'S', 'G3': 'F', 'G4': 'O', 'G5': 'P'}
        biom_add_metacol(obs, rankdic, 'Rank')
        exp = [{'Rank': 'S'}, {'Rank': 'S'}, {'Rank': 'F'}, {'Rank': 'O'},
               {'Rank': 'P'}]
        self.assertListEqual(list(map(
            dict, obs.metadata(axis='observation'))), exp)
        namedic = {'G1': 'Proteo', 'G3': 'Actino', 'G2': 'Firmic',
                   'G4': 'Bacter'}
        biom_add_metacol(obs, namedic, 'Name', missing='X')
        exp = [{'Rank': 'S', 'Name': 'Proteo'},
               {'Rank': 'S', 'Name': 'Firmic'},
               {'Rank': 'F', 'Name': 'Actino'},
               {'Rank': 'O', 'Name': 'Bacter'},
               {'Rank': 'P', 'Name': 'X'}]
        self.assertListEqual(list(map(
            dict, obs.metadata(axis='observation'))), exp)

    def test_collapse_biom(self):
        table = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8, 'G4': 0, 'G5': 3, 'G6': 0},
            'S2': {'G1': 1, 'G2': 8, 'G3': 0, 'G4': 7, 'G5': 4, 'G6': 2},
            'S3': {'G1': 0, 'G2': 2, 'G3': 3, 'G4': 5, 'G5': 0, 'G6': 9}})))

        # one-to-one mapping (e.g., direct translation)
        mapping = {'G1': ['H1'], 'G2': ['H2'], 'G3': ['H3'],
                   'G4': ['H4'], 'G5': ['H5'], 'G6': ['H6']}
        obs = collapse_biom(table.copy(), mapping)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 4, 'H2': 5, 'H3': 8, 'H4': 0, 'H5': 3, 'H6': 0},
            'S2': {'H1': 1, 'H2': 8, 'H3': 0, 'H4': 7, 'H5': 4, 'H6': 2},
            'S3': {'H1': 0, 'H2': 2, 'H3': 3, 'H4': 5, 'H5': 0, 'H6': 9}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # some missing, some extra
        mapping = {'G1': ['H1'], 'G2': ['H2'], 'G3': ['H3'], 'G9': ['H9']}
        obs = collapse_biom(table.copy(), mapping)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 4, 'H2': 5, 'H3': 8},
            'S2': {'H1': 1, 'H2': 8, 'H3': 0},
            'S3': {'H1': 0, 'H2': 2, 'H3': 3}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # wrong mapping (no match)
        mapping = {'H1': ['I1'], 'H2': ['I2'], 'H3': ['I3']}
        obs = collapse_biom(table.copy(), mapping)
        self.assertTrue(obs.is_empty())
        self.assertListEqual(list(obs.ids('sample')), ['S1', 'S2', 'S3'])
        self.assertListEqual(list(obs.ids('observation')), [])

        # many-to-one mapping (e.g., taxonomic rank up)
        mapping = {'G1': ['H1'], 'G2': ['H1'], 'G3': ['H2'],
                   'G4': ['H2'], 'G5': ['H2'], 'G6': ['H3']}
        obs = collapse_biom(table.copy(), mapping)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 9, 'H2': 11, 'H3': 0},
            'S2': {'H1': 9, 'H2': 11, 'H3': 2},
            'S3': {'H1': 2, 'H2':  8, 'H3': 9}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # many-to-many mapping (e.g., genes to pathways)
        mapping = {'G1': ['H1'],
                   'G2': ['H1', 'H2'],
                   'G3': ['H2', 'H3', 'H4'],
                   'G4': ['H2', 'H5'],
                   'G5': ['H4'],
                   'G6': ['H3', 'H5']}
        obs = collapse_biom(table.copy(), mapping)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 9, 'H2': 13, 'H3':  8, 'H4': 11, 'H5':  0},
            'S2': {'H1': 9, 'H2': 15, 'H3':  2, 'H4':  4, 'H5':  9},
            'S3': {'H1': 2, 'H2': 10, 'H3': 12, 'H4':  3, 'H5': 14}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # many-to-many mapping, with normalization
        obs = collapse_biom(table.copy(), mapping, normalize=True)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 6, 'H2': 5, 'H3': 3, 'H4': 6, 'H5': 0},
            'S2': {'H1': 5, 'H2': 8, 'H3': 1, 'H4': 4, 'H5': 4},
            'S3': {'H1': 1, 'H2': 4, 'H3': 6, 'H4': 1, 'H5': 7}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # nothing left after normalization
        table = Table(*map(np.array, prep_table({
            'S1': {'G1': 0}, 'S2': {'G1': 1}, 'S3': {'G1': 2}})))
        mapping = {'G1': ['H1', 'H2', 'H3', 'H4']}
        obs = collapse_biom(table.copy(), mapping, normalize=True)
        self.assertTrue(obs.is_empty())
        self.assertListEqual(list(obs.ids('sample')), ['S1', 'S2', 'S3'])
        self.assertListEqual(list(obs.ids('observation')), [])


if __name__ == '__main__':
    main()
