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
from biom import Table, load_table

from woltka.table import (
    prep_table, read_table, write_table, read_tsv, write_tsv, strip_metacols,
    table_shape, filter_table, merge_tables, round_table, add_metacol,
    collapse_table, calc_coverage)
from woltka.biom import table_to_biom


class TableTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_prep_table(self):
        # default mode
        prof = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                'S3': {'G2': 3, 'G5': 5}}
        obs = prep_table(prof)
        self.assertListEqual(obs[0], [
            [4, 2, 0], [5, 0, 3], [8, 0, 0], [0, 3, 0], [0, 7, 5]])
        self.assertListEqual(obs[1], ['G1', 'G2', 'G3', 'G4', 'G5'])
        self.assertListEqual(obs[2], ['S1', 'S2', 'S3'])
        self.assertListEqual(obs[3], [{}] * 5)

        # with sample Ids in custom order
        samples = ['S3', 'S1']
        obs = prep_table(prof, samples=samples)
        self.assertListEqual(obs[2], ['S3', 'S1'])
        self.assertListEqual(obs[0], [
            [0, 4], [3, 5], [0, 8], [5, 0]])

        # some sample Ids are not in data
        samples = ['S3', 'S0', 'S1']
        obs = prep_table(prof, samples=samples)
        self.assertListEqual(obs[2], ['S3', 'S1'])
        self.assertListEqual(obs[0], [
            [0, 4], [3, 5], [0, 8], [5, 0]])

        # with taxon names
        namedic = {'G1': 'Actinobacteria',
                   'G2': 'Firmicutes',
                   'G3': 'Bacteroidetes',
                   'G4': 'Cyanobacteria'}
        obs = prep_table(prof, namedic=namedic)
        self.assertListEqual(obs[1], ['G1', 'G2', 'G3', 'G4', 'G5'])
        self.assertListEqual([x['Name'] for x in obs[3]], [
            'Actinobacteria', 'Firmicutes', 'Bacteroidetes', 'Cyanobacteria',
            ''])

        # with taxon names to replace Ids
        obs = prep_table(prof, namedic=namedic, name_as_id=True)
        self.assertListEqual(obs[1], [
            'Actinobacteria', 'Firmicutes', 'Bacteroidetes', 'Cyanobacteria',
            'G5'])
        self.assertListEqual(obs[3], [{}] * 5)

        # with ranks
        rankdic = {'G1': 'class', 'G2': 'phylum', 'G4': 'phylum'}
        obs = prep_table(prof, rankdic=rankdic)
        self.assertListEqual([x['Rank'] for x in obs[3]], [
            'class', 'phylum', '', 'phylum', ''])

        # with lineages
        tree = {'G1': '74',  # Actinobacteria (phylum)
                '74': '72',
                'G2': '72',  # Terrabacteria group
                'G3': '70',  # FCB group
                'G4': '72',
                'G5': '1',
                '72': '2',
                '70': '2',
                '2':  '1',
                '1':  '1'}
        obs = prep_table(prof, tree=tree)
        self.assertListEqual([x['Lineage'] for x in obs[3]], [
            '2;72;74', '2;72', '2;70', '2;72', ''])

        # with lineages and names as Ids
        namedic.update({
            '74': 'Actino', '72': 'Terra', '70': 'FCB', '2': 'Bacteria'})
        obs = prep_table(prof, tree=tree, namedic=namedic, name_as_id=True)
        self.assertListEqual(obs[1], [
            'Actinobacteria', 'Firmicutes', 'Bacteroidetes', 'Cyanobacteria',
            'G5'])
        self.assertListEqual([x['Lineage'] for x in obs[3]], [
            'Bacteria;Terra;Actino', 'Bacteria;Terra', 'Bacteria;FCB',
            'Bacteria;Terra', ''])

        # with stratification
        sprof = {'S1': {('A', 'G1'): 4,
                        ('A', 'G2'): 5,
                        ('B', 'G1'): 8},
                 'S2': {('A', 'G1'): 2,
                        ('B', 'G1'): 3,
                        ('B', 'G2'): 7},
                 'S3': {('B', 'G3'): 3,
                        ('C', 'G2'): 5}}
        obs = prep_table(sprof)
        self.assertListEqual(obs[0], [
            [4, 2, 0], [5, 0, 0], [8, 3, 0], [0, 7, 0], [0, 0, 3], [0, 0, 5]])
        self.assertListEqual(obs[1], [
            'A|G1', 'A|G2', 'B|G1', 'B|G2', 'B|G3', 'C|G2'])
        self.assertListEqual(obs[2], ['S1', 'S2', 'S3'])

        # empty parameters instead of None
        obs = prep_table(prof, None, {}, {}, {})
        self.assertListEqual(obs[3], [{}] * 5)
        obs = prep_table(prof, [], {}, {}, {}, True)
        self.assertListEqual(obs[1], ['G1', 'G2', 'G3', 'G4', 'G5'])
        self.assertListEqual(obs[3], [{}] * 5)

    def test_read_table(self):
        # read a BIOM table
        fp = join(self.datdir, 'output', 'blastn.species.biom')
        table, fmt = read_table(fp)
        self.assertEqual(fmt, 'biom')

        # read a TSV file
        fp = join(self.datdir, 'output', 'blastn.species.tsv')
        table, fmt = read_table(fp)
        self.assertEqual(fmt, 'tsv')

        # wrong encoding
        fp = join(self.datdir, 'function', 'uniref', 'uniref.map.xz')
        with self.assertRaises(ValueError) as ctx:
            read_table(fp)
        errmsg = 'Input file cannot be parsed as BIOM or TSV format.'
        self.assertEqual(str(ctx.exception), errmsg)

        # error while parsing TSV
        fp = join(self.datdir, 'tree.nwk')
        with self.assertRaises(ValueError) as ctx:
            read_table(fp)
        errmsg = 'Input table file has no sample.'
        self.assertEqual(str(ctx.exception), errmsg)

    def test_write_table(self):
        table = (
            [[4, 2, 0],
             [5, 0, 3],
             [8, 0, 0],
             [0, 3, 0],
             [0, 7, 5]],
            ['G1', 'G2', 'G3', 'G4', 'G5'],
            ['S1', 'S2', 'S3'],
            [{'Name': 'Actinobacteria'},
             {'Name': 'Firmicutes'},
             {'Name': 'Bacteroidetes'},
             {'Name': 'Cyanobacteria'},
             {'Name': ''}])
        biota = table_to_biom(*table)

        # tuple to TSV
        fp = join(self.tmpdir, 'output.tsv')
        write_table(table, fp)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3\tName',
               'G1\t4\t2\t0\tActinobacteria',
               'G2\t5\t0\t3\tFirmicutes',
               'G3\t8\t0\t0\tBacteroidetes',
               'G4\t0\t3\t0\tCyanobacteria',
               'G5\t0\t7\t5\t']
        self.assertListEqual(obs, exp)

        # BIOM to TSV
        write_table(biota, fp)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertListEqual(obs, exp)
        remove(fp)

        # BIOM to BIOM
        fp = join(self.tmpdir, 'output.biom')
        write_table(biota, fp)
        obs = load_table(fp)
        self.assertEqual(obs.descriptive_equality(biota),
                         'Tables appear equal')

        # tuple to BIOM
        write_table(table, fp)
        obs = load_table(fp)
        self.assertEqual(obs.descriptive_equality(biota),
                         'Tables appear equal')
        remove(fp)

    def test_read_tsv(self):
        # data only
        tsv = ['#FeatureID\tS1\tS2\tS3',
               'G1\t4\t2\t0',
               'G2\t5\t0\t3',
               'G3\t8\t0\t0',
               'G4\t0\t3\t0',
               'G5\t0\t7\t5']
        obs = read_tsv(iter(tsv))
        self.assertListEqual(obs[0], [
            [4, 2, 0], [5, 0, 3], [8, 0, 0], [0, 3, 0], [0, 7, 5]])
        self.assertListEqual(obs[1], ['G1', 'G2', 'G3', 'G4', 'G5'])
        self.assertListEqual(obs[2], ['S1', 'S2', 'S3'])
        self.assertListEqual(obs[3], [{}] * 5)

        # with metadata
        tsv = ['#FeatureID\tS1\tS2\tS3\tName\tRank\tLineage',
               'G1\t4\t2\t0\tActinobacteria\tphylum\t2;72;74',
               'G2\t5\t0\t3\tFirmicutes\tphylum\t2;72',
               'G3\t8\t0\t0\tBacteroidetes\tphylum\t2;70',
               'G4\t0\t3\t0\tCyanobacteria\tphylum\t2;72',
               'G5\t0\t7\t5\t\t\t']
        obs = read_tsv(iter(tsv))
        self.assertListEqual(obs[0], [
            [4, 2, 0], [5, 0, 3], [8, 0, 0], [0, 3, 0], [0, 7, 5]])
        self.assertListEqual(obs[1], ['G1', 'G2', 'G3', 'G4', 'G5'])
        self.assertListEqual(obs[2], ['S1', 'S2', 'S3'])
        self.assertListEqual(obs[3], [
            {'Name': 'Actinobacteria', 'Rank': 'phylum', 'Lineage': '2;72;74'},
            {'Name': 'Firmicutes',     'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': 'Bacteroidetes',  'Rank': 'phylum', 'Lineage': '2;70'},
            {'Name': 'Cyanobacteria',  'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': '',               'Rank': '',       'Lineage': ''}])

        # empty file
        with self.assertRaises(ValueError) as ctx:
            read_tsv(iter([]))
        self.assertEqual(str(ctx.exception), 'Input table file is empty.')

        # no sample
        with self.assertRaises(ValueError) as ctx:
            read_tsv(iter(['#ID\tName']))
        self.assertEqual(str(ctx.exception), 'Input table file has no sample.')

    def test_write_tsv(self):
        fp = join(self.tmpdir, 'table.tsv')

        # just data
        data = [[4, 2, 0],
                [5, 0, 3],
                [8, 0, 0],
                [0, 3, 0],
                [0, 7, 5]]
        features = ['G1', 'G2', 'G3', 'G4', 'G5']
        samples = ['S1', 'S2', 'S3']
        with open(fp, 'w') as f:
            write_tsv((data, features, samples, None), f)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3',
               'G1\t4\t2\t0',
               'G2\t5\t0\t3',
               'G3\t8\t0\t0',
               'G4\t0\t3\t0',
               'G5\t0\t7\t5']
        self.assertListEqual(obs, exp)

        # with metadata
        metadata = [
            {'Name': 'Actinobacteria', 'Rank': 'phylum', 'Lineage': '2;72;74'},
            {'Name': 'Firmicutes',     'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': 'Bacteroidetes',  'Rank': 'phylum', 'Lineage': '2;70'},
            {'Name': 'Cyanobacteria',  'Rank': 'phylum', 'Lineage': '2;72'},
            {'Name': '',               'Rank': '',       'Lineage': ''}]
        with open(fp, 'w') as f:
            write_tsv((data, features, samples, metadata), f)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = [
            '#FeatureID\tS1\tS2\tS3\tName\tRank\tLineage',
            'G1\t4\t2\t0\tActinobacteria\tphylum\t2;72;74',
            'G2\t5\t0\t3\tFirmicutes\tphylum\t2;72',
            'G3\t8\t0\t0\tBacteroidetes\tphylum\t2;70',
            'G4\t0\t3\t0\tCyanobacteria\tphylum\t2;72',
            'G5\t0\t7\t5\t\t\t']
        self.assertListEqual(obs, exp)
        remove(fp)

    def test_strip_cols(self):
        # all three metadata columns
        header = ['#ID', 'S01', 'S02', 'S03', 'Name', 'Rank', 'Lineage']
        obs = strip_metacols(header)
        self.assertListEqual(obs[0], ['#ID', 'S01', 'S02', 'S03'])
        self.assertListEqual(obs[1], ['Name', 'Rank', 'Lineage'])

        # no metadata column
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'S01', 'S03']),
            (['#ID', 'S01', 'S01', 'S03'], []))

        # 1st column
        self.assertTupleEqual(strip_metacols(['#ID', 'S01', 'Name']), (
            ['#ID', 'S01'], ['Name']))

        # 2nd column
        self.assertTupleEqual(strip_metacols(['#ID', 'S01', 'Rank']), (
            ['#ID', 'S01'], ['Rank']))

        # last column
        self.assertTupleEqual(strip_metacols(['#ID', 'S01', 'Lineage']), (
            ['#ID', 'S01'], ['Lineage']))

        # 1st and 2nd columns
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'Name', 'Rank']),
            (['#ID', 'S01'], ['Name', 'Rank']))

        # 1st and last columns
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'Name', 'Lineage']),
            (['#ID', 'S01'], ['Name', 'Lineage']))

        # 2nd and last columns
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'Rank', 'Lineage']),
            (['#ID', 'S01'], ['Rank', 'Lineage']))

        # only metadata columns
        self.assertTupleEqual(strip_metacols(['Name', 'Rank', 'Lineage']), (
            [], ['Name', 'Rank', 'Lineage']))

        # metadata column mixed in samples (ignored)
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'S01', 'Name', 'S03']),
            (['#ID', 'S01', 'S01', 'Name', 'S03'], []))

        # metadata column in wrong order
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'Lineage', 'Name', 'Rank']),
            (['#ID', 'S01', 'Lineage'], ['Name', 'Rank']))

        # duplicate metadata column
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'Name', 'Name', 'Rank']),
            (['#ID', 'S01', 'Name'], ['Name', 'Rank']))

        # duplicate metadata column in wrong order
        self.assertTupleEqual(
            strip_metacols(['#ID', 'S01', 'Rank', 'Lineage', 'Rank']),
            (['#ID', 'S01', 'Rank', 'Lineage'], ['Rank']))

        # nothing
        self.assertTupleEqual(strip_metacols([]), ([], []))
        self.assertTupleEqual(strip_metacols([], []), ([], []))

        # custom metadata columns
        header = ['#ID', 'S01', 'S02', 'S03', 'Name', 'Rank', 'Lineage']
        obs = strip_metacols(header, ['Rank', 'Lineage'])
        self.assertListEqual(obs[0], ['#ID', 'S01', 'S02', 'S03', 'Name'])
        self.assertListEqual(obs[1], ['Rank', 'Lineage'])

    def test_table_shape(self):
        # original tuple
        table = prep_table({'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                            'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                            'S3': {'G2': 3, 'G5': 5}})
        self.assertTupleEqual(table_shape(table), (5, 3))

        # BIOM table
        table = Table(*map(np.array, table))
        self.assertTupleEqual(table_shape(table), (5, 3))

    def test_filter_table(self):
        table = prep_table({'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                            'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                            'S3': {'G2': 3, 'G5': 5}})

        # filter by count
        obs = filter_table(table, th=3)
        exp = ([[4, 0, 0], [5, 0, 3], [8, 0, 0], [0, 3, 0], [0, 7, 5]],
               ['G1', 'G2', 'G3', 'G4', 'G5'], ['S1', 'S2', 'S3'], [{}] * 5)
        self.assertTupleEqual(obs, exp)

        obs = filter_table(table, th=4)
        exp = ([[4, 0, 0], [5, 0, 0], [8, 0, 0], [0, 7, 5]],
               ['G1', 'G2', 'G3', 'G5'], ['S1', 'S2', 'S3'], [{}] * 4)
        self.assertTupleEqual(obs, exp)

        obs = filter_table(table, th=6)
        exp = ([[8, 0, 0], [0, 7, 0]], ['G3', 'G5'], ['S1', 'S2', 'S3'],
               [{}] * 2)
        self.assertTupleEqual(obs, exp)

        # filter by threshold
        obs = filter_table(table, th=0.25)
        exp = ([[5, 0, 3], [8, 0, 0], [0, 3, 0], [0, 7, 5]],
               ['G2', 'G3', 'G4', 'G5'], ['S1', 'S2', 'S3'], [{}] * 4)
        self.assertTupleEqual(obs, exp)

        obs = filter_table(table, th=0.5)
        exp = ([[0, 7, 5]], ['G5'], ['S1', 'S2', 'S3'], [{}])
        self.assertTupleEqual(obs, exp)

        # filter out everything
        obs = filter_table(table, th=10)
        exp = ([], [], ['S1', 'S2', 'S3'], [])
        self.assertTupleEqual(obs, exp)

        # filter an empty table
        obs = filter_table(exp, th=1)
        exp = ([], [], ['S1', 'S2', 'S3'], [])
        self.assertTupleEqual(obs, exp)

        # filter a BIOM table
        table = Table(*map(np.array, table))
        obs = filter_table(table, th=3)
        exp = Table(*map(np.array, prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8},
            'S2': {'G4': 3, 'G5': 7},
            'S3': {'G2': 3, 'G5': 5}})))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

    def test_merge_tables(self):
        # just data
        t1 = prep_table({'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                         'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                         'S3': {'G2': 3, 'G5': 5}})
        t2 = prep_table({'S3': {'G3': 1, 'G5': 1},
                         'S4': {'G2': 5, 'G3': 3, 'G6': 9},
                         'S5': {'G5': 2, 'G6': 4}})
        t3 = prep_table({'S2': {'G3': 2, 'G5': 2, 'G6': 8},
                         'S6': {'G3': 1, 'G6': 6}})
        obs = merge_tables([t1, t2, t3])
        exp = prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8, 'G4': 0, 'G5': 0, 'G6': 0},
            'S2': {'G1': 2, 'G2': 0, 'G3': 2, 'G4': 3, 'G5': 9, 'G6': 8},
            'S3': {'G1': 0, 'G2': 3, 'G3': 1, 'G4': 0, 'G5': 6, 'G6': 0},
            'S4': {'G1': 0, 'G2': 5, 'G3': 3, 'G4': 0, 'G5': 0, 'G6': 9},
            'S5': {'G1': 0, 'G2': 0, 'G3': 0, 'G4': 0, 'G5': 2, 'G6': 4},
            'S6': {'G1': 0, 'G2': 0, 'G3': 1, 'G4': 0, 'G5': 0, 'G6': 6}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # with metadata
        names = {'G1': 'Actinobacteria',
                 'G2': 'Firmicutes',
                 'G3': 'Bacteroidetes',
                 'G4': 'Cyanobacteria',
                 'G5': 'Proteobacteria',
                 'G6': 'Fusobacteria'}
        for t in (t1, t2, t3, exp):
            t[3].clear()
            t[3].extend({'Name': names[x]} for x in t[1])
        obs = merge_tables([t1, t2, t3])
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # some biom tables
        obs = merge_tables([t1, table_to_biom(*t2), t3])
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # all biom tables
        obs = merge_tables([table_to_biom(*x) for x in (t1, t2, t3)])
        self.assertTrue(isinstance(obs, Table))
        exp = table_to_biom(*exp)
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # inconsistent metadata
        t3[3][1]['Name'] = 'This is not right.'
        with self.assertRaises(ValueError) as ctx:
            merge_tables([t1, t2, t3])
        errmsg = 'Conflicting metadata found in tables.'
        self.assertEqual(str(ctx.exception), errmsg)

    def test_round_table(self):
        obs = prep_table({
            'S1': {'G1': 0.5, 'G2': 0.0, 'G3': 2.3, 'G4': 0.50000000001},
            'S2': {'G1': 1.5, 'G2': 0.2, 'G3': 1.49999999999, 'G4': 0.2},
            'S3': {'G1': 2.5, 'G2': 0.3, 'G3': 3.8, 'G4': 0.1}})
        exp = prep_table({
            'S1': {'G1': 0, 'G3': 2},
            'S2': {'G1': 2, 'G3': 2},
            'S3': {'G1': 2, 'G3': 4}})
        ob2 = Table(*map(np.array, obs))

        # regular
        round_table(obs)
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # BIOM
        round_table(ob2)
        ex2 = Table(*map(np.array, exp))
        self.assertEqual(ob2.descriptive_equality(ex2), 'Tables appear equal')

    def test_add_metacol(self):
        obs = prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8, 'G4': 0, 'G5': 3},
            'S2': {'G1': 1, 'G2': 8, 'G3': 0, 'G4': 7, 'G5': 4},
            'S3': {'G1': 0, 'G2': 2, 'G3': 3, 'G4': 5, 'G5': 0}})
        self.assertListEqual(obs[3], [{}] * 5)
        ob2 = Table(*map(np.array, obs))

        # regular table
        rankdic = {'G1': 'S', 'G2': 'S', 'G3': 'F', 'G4': 'O', 'G5': 'P'}
        add_metacol(obs, rankdic, 'Rank')
        exp = [{'Rank': 'S'}, {'Rank': 'S'}, {'Rank': 'F'}, {'Rank': 'O'},
               {'Rank': 'P'}]
        self.assertListEqual(obs[3], exp)

        # BIOM table
        add_metacol(ob2, rankdic, 'Rank')
        self.assertListEqual(list(map(
            dict, ob2.metadata(axis='observation'))), exp)

        # unordered, missing value, append
        namedic = {'G1': 'Proteo', 'G3': 'Actino', 'G2': 'Firmic',
                   'G4': 'Bacter'}
        add_metacol(obs, namedic, 'Name', missing='X')
        exp = [{'Rank': 'S', 'Name': 'Proteo'},
               {'Rank': 'S', 'Name': 'Firmic'},
               {'Rank': 'F', 'Name': 'Actino'},
               {'Rank': 'O', 'Name': 'Bacter'},
               {'Rank': 'P', 'Name': 'X'}]
        self.assertListEqual(obs[3], exp)

        add_metacol(ob2, namedic, 'Name', missing='X')
        self.assertListEqual(list(map(
            dict, ob2.metadata(axis='observation'))), exp)

    def test_collapse_table(self):
        table = prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8, 'G4': 0, 'G5': 3, 'G6': 0},
            'S2': {'G1': 1, 'G2': 8, 'G3': 0, 'G4': 7, 'G5': 4, 'G6': 2},
            'S3': {'G1': 0, 'G2': 2, 'G3': 3, 'G4': 5, 'G5': 0, 'G6': 9}})

        # one-to-one mapping (e.g., direct translation)
        mapping = {'G1': ['H1'], 'G2': ['H2'], 'G3': ['H3'],
                   'G4': ['H4'], 'G5': ['H5'], 'G6': ['H6']}
        obs = collapse_table(table, mapping)
        exp = prep_table({
            'S1': {'H1': 4, 'H2': 5, 'H3': 8, 'H4': 0, 'H5': 3, 'H6': 0},
            'S2': {'H1': 1, 'H2': 8, 'H3': 0, 'H4': 7, 'H5': 4, 'H6': 2},
            'S3': {'H1': 0, 'H2': 2, 'H3': 3, 'H4': 5, 'H5': 0, 'H6': 9}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # BIOM table
        table_ = Table(*map(np.array, table))
        obs = collapse_table(table_, mapping)
        exp = Table(*map(np.array, exp))
        self.assertEqual(obs.descriptive_equality(exp), 'Tables appear equal')

        # some missing, some extra
        mapping = {'G1': ['H1'], 'G2': ['H2'], 'G3': ['H3'], 'G9': ['H9']}
        obs = collapse_table(table, mapping)
        exp = prep_table({
            'S1': {'H1': 4, 'H2': 5, 'H3': 8},
            'S2': {'H1': 1, 'H2': 8, 'H3': 0},
            'S3': {'H1': 0, 'H2': 2, 'H3': 3}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # wrong mapping (no match)
        mapping = {'H1': ['I1'], 'H2': ['I2'], 'H3': ['I3']}
        obs = collapse_table(table, mapping)
        for i in (0, 1, 3):
            self.assertListEqual(obs[i], [])
        self.assertListEqual(obs[2], ['S1', 'S2', 'S3'])

        # many-to-one mapping (e.g., taxonomic rank up)
        mapping = {'G1': ['H1'], 'G2': ['H1'], 'G3': ['H2'],
                   'G4': ['H2'], 'G5': ['H2'], 'G6': ['H3']}
        obs = collapse_table(table, mapping)
        exp = prep_table({
            'S1': {'H1': 9, 'H2': 11, 'H3': 0},
            'S2': {'H1': 9, 'H2': 11, 'H3': 2},
            'S3': {'H1': 2, 'H2':  8, 'H3': 9}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # many-to-many mapping (e.g., genes to pathways)
        mapping = {'G1': ['H1'],
                   'G2': ['H1', 'H2'],
                   'G3': ['H2', 'H3', 'H4'],
                   'G4': ['H2', 'H5'],
                   'G5': ['H4'],
                   'G6': ['H3', 'H5']}
        obs = collapse_table(table, mapping)
        exp = prep_table({
            'S1': {'H1': 9, 'H2': 13, 'H3':  8, 'H4': 11, 'H5':  0},
            'S2': {'H1': 9, 'H2': 15, 'H3':  2, 'H4':  4, 'H5':  9},
            'S3': {'H1': 2, 'H2': 10, 'H3': 12, 'H4':  3, 'H5': 14}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # many-to-many mapping, with normalization
        obs = collapse_table(table, mapping, normalize=True)
        exp = prep_table({
            'S1': {'H1': 6, 'H2': 5, 'H3': 3, 'H4': 6, 'H5': 0},
            'S2': {'H1': 5, 'H2': 8, 'H3': 1, 'H4': 4, 'H5': 4},
            'S3': {'H1': 1, 'H2': 4, 'H3': 6, 'H4': 1, 'H5': 7}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # nothing left after normalization
        table = prep_table({'S1': {'G1': 0}, 'S2': {'G1': 1}, 'S3': {'G1': 2}})
        mapping = {'G1': ['H1', 'H2', 'H3', 'H4']}
        obs = collapse_table(table, mapping, normalize=True)
        for i in (0, 1, 3):
            self.assertListEqual(obs[i], [])
        self.assertListEqual(obs[2], ['S1', 'S2', 'S3'])

    def test_calc_coverage(self):
        table = prep_table({
            'S1': {'G1': 4, 'G2': 5, 'G3': 8, 'G4': 0, 'G5': 3, 'G6': 0},
            'S2': {'G1': 1, 'G2': 8, 'G3': 0, 'G4': 7, 'G5': 4, 'G6': 2},
            'S3': {'G1': 0, 'G2': 2, 'G3': 3, 'G4': 5, 'G5': 0, 'G6': 9}})
        mapping = {'P1': ['G1', 'G2'],
                   'P2': ['G3'],
                   'P3': ['G2', 'G4', 'G6'],
                   'P4': ['G3', 'G5'],
                   'P5': ['G7', 'G8', 'G9']}

        # default behavior
        obs = calc_coverage(table, mapping)
        exp = prep_table({
            'S1': {'P1': 100.0, 'P2': 100.0, 'P3':  33.333, 'P4': 100.0},
            'S2': {'P1': 100.0, 'P2':   0.0, 'P3': 100.0,   'P4': 50.0},
            'S3': {'P1':  50.0, 'P2': 100.0, 'P3': 100.0,   'P4': 50.0}})
        for i in range(4):
            self.assertListEqual(obs[i], exp[i])

        # BIOM table
        table_ = Table(*map(np.array, table))
        obs = calc_coverage(table_, mapping)
        for i in range(2):
            self.assertListEqual(obs[i], exp[i])

        # threshold and boolean result
        obs = calc_coverage(table, mapping, th=80)
        exp = prep_table({
            'S1': {'P1': 1, 'P2': 1, 'P3': 0, 'P4': 1},
            'S2': {'P1': 1, 'P2': 0, 'P3': 1, 'P4': 0},
            'S3': {'P1': 0, 'P2': 1, 'P3': 1, 'P4': 0}})
        for i in range(2):
            self.assertListEqual(obs[i], exp[i])

        # numbers instead of percentages
        obs = calc_coverage(table, mapping, count=True)
        exp = prep_table({
            'S1': {'P1': 2, 'P2': 1, 'P3': 1, 'P4': 2},
            'S2': {'P1': 2, 'P2': 0, 'P3': 3, 'P4': 1},
            'S3': {'P1': 1, 'P2': 1, 'P3': 3, 'P4': 1}})
        for i in range(2):
            self.assertListEqual(obs[i], exp[i])

        # number overrides threshold
        obs = calc_coverage(table, mapping, th=80, count=True)
        for i in range(2):
            self.assertListEqual(obs[i], exp[i])


if __name__ == '__main__':
    main()
