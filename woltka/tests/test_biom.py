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
    table_to_biom, biom_to_table, write_biom, biom_max_f, divide_biom,
    scale_biom, filter_biom, round_biom, biom_add_metacol, clip_biom,
    collapse_biom)
from woltka.table import prep_table


class BiomTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def assertBIOMEqual(self, a, b):
        """Assert two BIOM tables are equal.
        """
        self.assertEqual(a.descriptive_equality(b), 'Tables appear equal')

    def assertBIOMEmpty(self, a, samples=None):
        """Assert a BIOM table is empty; optionally check sample IDs.
        """
        self.assertTrue(a.is_empty())
        self.assertListEqual(list(a.ids('observation')), [])
        if samples:
            self.assertListEqual(list(a.ids('sample')), samples)

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
        self.assertBIOMEqual(obs, exp)
        remove(fp)

    def test_biom_max_f(self):
        table = Table(*map(np.array, prep_table({
            'S1': {'G1': 0, 'G2': 1, 'G3': 200},
            'S2': {'G1': 1.5, 'G2': 2.475, 'G3': 8.12782},
            'S3': {'G1': 1e-5, 'G2': 33.905, 'G3': 3.1415926}})))
        self.assertEqual(biom_max_f(table), 7)

    def test_divide_biom(self):
        obs = Table(*map(np.array, prep_table({
            'S1': {'G1': 20, 'G2': 36, 'G3': 4},
            'S2': {'G1': 15, 'G2': 24, 'G3': 8},
            'S3': {'G1': 10, 'G2': 18, 'G3': 0}})))
        sizes = {'G1': 5, 'G2': 6, 'G3': 2}
        divide_biom(obs, sizes)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 6, 'G3': 2},
            'S2': {'G1': 3, 'G2': 4, 'G3': 4},
            'S3': {'G1': 2, 'G2': 3, 'G3': 0}})))
        self.assertBIOMEqual(obs, exp)
        del sizes['G3']
        with self.assertRaises(KeyError):
            divide_biom(obs, sizes)

    def test_scale_biom(self):
        obs = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 7, 'G3': 0},
            'S2': {'G1': 2, 'G2': 3, 'G3': 1}})))
        scale_biom(obs, 3)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 12, 'G2': 21, 'G3': 0},
            'S2': {'G1':  6, 'G2':  9, 'G3': 3}})))
        self.assertBIOMEqual(obs, exp)

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
        self.assertBIOMEqual(obs, exp)

        obs = filter_biom(table, th=4)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8},
            'S2': {'G5': 7},
            'S3': {'G5': 5}})))
        self.assertBIOMEqual(obs, exp)

        obs = filter_biom(table, th=6)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G3': 8},
            'S2': {'G5': 7},
            'S3': {}})))
        self.assertBIOMEqual(obs, exp)

        obs = filter_biom(table, th=0.25)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G2': 5, 'G3': 8},
            'S2': {'G4': 3, 'G5': 7},
            'S3': {'G2': 3, 'G5': 5}})))
        self.assertBIOMEqual(obs, exp)

        obs = filter_biom(table, th=0.5)
        exp = Table(*map(np.array, prep_table({
            'S1': {},
            'S2': {'G5': 7},
            'S3': {'G5': 5}})))
        self.assertBIOMEqual(obs, exp)

        obs = filter_biom(table, th=10)
        self.assertBIOMEmpty(obs, ['S1', 'S2', 'S3'])

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
        self.assertBIOMEqual(obs, exp)
        obs = Table(*map(np.array, prep_table({
            'S1': {'G1': 0.225, 'G2': 0.0,   'G3': 2.375},
            'S2': {'G1': 1.547, 'G2': 0.173, 'G3': 1.499}})))
        round_biom(obs, 2)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 0.23, 'G2': 0.0,  'G3': 2.38},
            'S2': {'G1': 1.55, 'G2': 0.17, 'G3': 1.5}})))
        self.assertBIOMEqual(obs, exp)

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

    def test_clip_biom(self):
        table = Table(*map(np.array, prep_table({
            'S1': {'G1_1': 4, 'G1_2': 5, 'G1_3': 0, 'G2_1': 0, 'G2_2': 3},
            'S2': {'G1_1': 1, 'G1_2': 8, 'G1_4': 0, 'G2_1': 3, 'G2_3': 4},
            'S3': {'G1_1': 0, 'G1_3': 2, 'G1_4': 3, 'G2_2': 5, 'G2_3': 0}})))

        # 1st field
        obs = clip_biom(table.copy(), 1, sep='_')
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 9, 'G2': 3},
            'S2': {'G1': 9, 'G2': 7},
            'S3': {'G1': 5, 'G2': 5}})))
        self.assertBIOMEqual(obs, exp)

        # nested
        obs = clip_biom(table.copy(), 1, sep='_', nested=True)
        self.assertBIOMEqual(obs, exp)

        # 2nd field
        obs = clip_biom(table.copy(), 2, sep='_')
        exp = Table(*map(np.array, prep_table({
            'S1': {'1': 4, '2': 8, '3': 0, '4': 0},
            'S2': {'1': 4, '2': 8, '3': 4, '4': 0},
            'S3': {'1': 0, '2': 5, '3': 2, '4': 3}})))
        self.assertBIOMEqual(obs, exp)

        # nested (no change)
        obs = clip_biom(table.copy(), 2, sep='_', nested=True)
        self.assertBIOMEqual(obs, table)

        # field number too large (all dropped)
        obs = clip_biom(table.copy(), 3, sep='_')
        self.assertBIOMEmpty(obs, ['S1', 'S2', 'S3'])

        # invalid separator
        obs = clip_biom(table.copy(), 1, sep='.')
        self.assertBIOMEqual(obs, table)

        # empty fields
        table = Table(*map(np.array, prep_table({
            'S1': {'_G1_1': 3, 'G2__3': 5, 'G5_4_': 1, '__G0': 2}})))
        obs = clip_biom(table.copy(), 2, sep='_')
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 3, '4': 1}})))
        self.assertBIOMEqual(obs, exp)

        # nested
        obs = clip_biom(table.copy(), 2, sep='_', nested=True)
        exp = Table(*map(np.array, prep_table({
            'S1': {'_G1': 3, 'G5_4': 1}})))
        self.assertBIOMEqual(obs, exp)

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
        self.assertBIOMEqual(obs, exp)

        # some missing, some extra
        mapping = {'G1': ['H1'], 'G2': ['H2'], 'G3': ['H3'], 'G9': ['H9']}
        obs = collapse_biom(table.copy(), mapping)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 4, 'H2': 5, 'H3': 8},
            'S2': {'H1': 1, 'H2': 8, 'H3': 0},
            'S3': {'H1': 0, 'H2': 2, 'H3': 3}})))
        self.assertBIOMEqual(obs, exp)

        # wrong mapping (no match)
        mapping = {'H1': ['I1'], 'H2': ['I2'], 'H3': ['I3']}
        obs = collapse_biom(table.copy(), mapping)
        self.assertBIOMEmpty(obs, ['S1', 'S2', 'S3'])

        # many-to-one mapping (e.g., taxonomic rank up)
        mapping = {'G1': ['H1'], 'G2': ['H1'], 'G3': ['H2'],
                   'G4': ['H2'], 'G5': ['H2'], 'G6': ['H3']}
        obs = collapse_biom(table.copy(), mapping)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 9, 'H2': 11, 'H3': 0},
            'S2': {'H1': 9, 'H2': 11, 'H3': 2},
            'S3': {'H1': 2, 'H2':  8, 'H3': 9}})))
        self.assertBIOMEqual(obs, exp)

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
        self.assertBIOMEqual(obs, exp)

        # many-to-many mapping, with normalization
        obs = collapse_biom(table.copy(), mapping, divide=True)
        exp = Table(*map(np.array, prep_table({
            'S1': {'H1': 6, 'H2': 5, 'H3': 3, 'H4': 6, 'H5': 0},
            'S2': {'H1': 5, 'H2': 8, 'H3': 1, 'H4': 4, 'H5': 4},
            'S3': {'H1': 1, 'H2': 4, 'H3': 6, 'H4': 1, 'H5': 7}})))
        self.assertBIOMEqual(obs, exp)

        # nothing left after normalization
        table = Table(*map(np.array, prep_table({
            'S1': {'G1': 0}, 'S2': {'G1': 1}, 'S3': {'G1': 2}})))
        mapping = {'G1': ['H1', 'H2', 'H3', 'H4']}
        obs = collapse_biom(table.copy(), mapping, divide=True)
        self.assertBIOMEmpty(obs, ['S1', 'S2', 'S3'])

        # stratified features
        table = Table(*map(np.array, prep_table({
            'S1': {'A|K1': 4, 'A|K2': 5, 'B|K2': 8, 'C|K3': 3, 'C|K4': 0},
            'S2': {'A|K1': 1, 'A|K2': 8, 'B|K2': 0, 'C|K3': 4, 'C|K4': 2}})))
        mapping = {'A': ['1'], 'B': ['1']}
        obs = collapse_biom(table.copy(), mapping, field=1, sep='|')
        exp = Table(*map(np.array, prep_table({
            'S1': {'1|K1': 4, '1|K2': 13},
            'S2': {'1|K1': 1, '1|K2': 8}})))
        self.assertBIOMEqual(obs, exp)
        mapping = {'K1': ['H1'], 'K2': ['H2', 'H3'], 'K3': ['H3']}
        obs = collapse_biom(table.copy(), mapping, field=2, sep='|')
        exp = Table(*map(np.array, prep_table({
            'S1': {'A|H1': 4, 'A|H2': 5, 'A|H3': 5, 'B|H2': 8, 'B|H3': 8,
                   'C|H3': 3},
            'S2': {'A|H1': 1, 'A|H2': 8, 'A|H3': 8, 'B|H2': 0, 'B|H3': 0,
                   'C|H3': 4}})))
        self.assertBIOMEqual(obs, exp)

        # invalid or empty field
        table = Table(*map(np.array, prep_table({
            'S1': {'G_1': 6, '||G2': 3},
            'S2': {'G|1': 1, 'G2|': 7}})))
        mapping = {'G1': ['H1'], 'G2': ['H2']}
        obs = collapse_biom(table.copy(), mapping, field=2, sep='|')
        self.assertBIOMEmpty(obs, ['S1', 'S2'])

        # nested features - 1st level
        table = Table(*map(np.array, prep_table({
            'S1': {'A_1': 3, 'A_2': 6, 'B_1': 7, 'B_2': 0},
            'S2': {'A_2': 2, 'B_3': 2, 'C_1': 4, 'C_3': 2}})))
        mapping = {'A': ['X'], 'B': ['X'], 'C': ['Y']}
        obs = collapse_biom(table.copy(), mapping, field=1, sep='_',
                            nested=True)
        exp = Table(*map(np.array, prep_table({
            'S1': {'X|A_1': 3, 'X|A_2': 6, 'X|B_1': 7, 'Y|B_2': 0},
            'S2': {'X|A_2': 2, 'X|B_3': 2, 'Y|C_1': 4, 'Y|C_3': 2}})))
        self.assertBIOMEqual(obs, exp)

        # 2nd level
        mapping = {'A_1': ['a'], 'A_2': ['b'],
                   'B_1': ['a'], 'B_2': ['b'],
                   'C_1': ['a'], 'C_2': ['b']}
        obs = collapse_biom(table.copy(), mapping, field=2, sep='_',
                            nested=True)
        exp = Table(*map(np.array, prep_table({
            'S1': {'A|a': 3, 'A|b': 6, 'B|a': 7, 'B|b': 0},
            'S2': {'A|b': 2, 'C|a': 4}})))
        self.assertBIOMEqual(obs, exp)


if __name__ == '__main__':
    main()
