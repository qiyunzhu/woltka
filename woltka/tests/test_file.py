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
import gzip
import bz2

from woltka.file import (
    openzip, file2stem, path2stem, read_ids, id2file_from_dir,
    id2file_from_map, read_map, write_readmap, write_table, prep_table)


class FileTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_openzip(self):
        text = 'Hello World!'

        # read regular file
        fp = join(self.tmpdir, 'test.txt')
        with open(fp, 'w') as f:
            f.write(text)
        with openzip(fp) as f:
            obs = f.read()
        self.assertEqual(obs, text)

        # write regular file
        with openzip(fp, 'wt') as f:
            f.write('Here I am!')
        with open(fp, 'r') as f:
            obs = f.read()
        self.assertEqual(obs, 'Here I am!')
        remove(fp)

        # read compressed file
        fpz = join(self.tmpdir, 'test.txt.gz')
        with gzip.open(fpz, 'wb') as f:
            f.write(text.encode())
        with openzip(fpz) as f:
            obs = f.read()
        self.assertEqual(obs, text)
        remove(fpz)

        # write compressed file
        fpb = join(self.tmpdir, 'test.txt.bz2')
        with openzip(fpb, 'at') as f:
            f.write('Here I am!')
        with bz2.open(fpb, 'rt') as f:
            obs = f.read()
        self.assertEqual(obs, 'Here I am!')
        remove(fpb)

    def test_file2stem(self):
        self.assertEqual(file2stem('input.txt'), 'input')
        self.assertEqual(file2stem('input.gz'), 'input')
        self.assertEqual(file2stem('input.txt.gz'), 'input')
        self.assertEqual(file2stem('input.txt.gz', ext='.gz'), 'input.txt')
        with self.assertRaises(ValueError) as ctx:
            file2stem('input.txt', ext='.gz')
        self.assertEqual(str(ctx.exception), (
            'Filepath and filename extension do not match.'))

    def test_path2stem(self):
        self.assertEqual(path2stem('input.txt.bz2'), 'input')
        self.assertEqual(path2stem('/home/input.txt'), 'input')
        self.assertEqual(path2stem('/home/input.fa', ext='a'), 'input.f')
        self.assertEqual(path2stem('/home/.bashrc'), '.bashrc')

    def test_read_ids(self):
        # simple list
        ids = ['a', 'b', 'c']
        obs = read_ids(ids)
        self.assertListEqual(obs, ids)

        # metadata table
        lines = ['#SampleID\tHeight\tWeight',
                 'S01\t120.5\t56.8',
                 'S02\t114.2\t52.5',
                 'S03\t117.8\t55.0',
                 '']
        obs = read_ids(lines)
        self.assertListEqual(obs, ['S01', 'S02', 'S03'])

        # no Id
        with self.assertRaises(ValueError) as ctx:
            read_ids(['#'])
        self.assertEqual(str(ctx.exception), 'No ID is read.')

        # duplicate Ids
        with self.assertRaises(ValueError) as ctx:
            read_ids(['S01', 'S02', 'S01', 'S03'])
        self.assertEqual(str(ctx.exception), 'Duplicate IDs found.')

        # nothing
        self.assertIsNone(read_ids(None))

    def test_id2file_from_dir(self):
        ids = ['a', 'b', 'c']

        # regular Fasta files
        for id_ in ids:
            open(join(self.tmpdir, '{}.faa'.format(id_)), 'a').close()
        obs = id2file_from_dir(self.tmpdir)
        exp = {x: '{}.faa'.format(x) for x in ids}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, '{}.faa'.format(id_)))

        # gzipped Fasta files
        for id_ in ids:
            open(join(self.tmpdir, '{}.faa.gz'.format(id_)), 'a').close()
        obs = id2file_from_dir(self.tmpdir)
        exp = {x: '{}.faa.gz'.format(x) for x in ids}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, '{}.faa.gz'.format(id_)))

        # user-defined extension filename
        for id_ in ids:
            open(join(self.tmpdir, '{}.faa'.format(id_)), 'a').close()
        open(join(self.tmpdir, 'readme.txt'), 'a').close()
        open(join(self.tmpdir, 'taxdump'), 'a').close()
        obs = id2file_from_dir(self.tmpdir, ext='.faa')
        exp = {x: '{}.faa'.format(x) for x in ids}
        self.assertDictEqual(obs, exp)
        remove(join(self.tmpdir, 'readme.txt'))
        remove(join(self.tmpdir, 'taxdump'))

        # user-defined Id list
        obs = id2file_from_dir(self.tmpdir, ids=['a', 'b'])
        exp = {x: '{}.faa'.format(x) for x in ['a', 'b']}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, '{}.faa'.format(id_)))

        # duplicated Ids
        open(join(self.tmpdir, 'x.faa'), 'a').close()
        open(join(self.tmpdir, 'x.tsv'), 'a').close()
        with self.assertRaises(ValueError) as ctx:
            id2file_from_dir(self.tmpdir)
        msg = 'Ambiguous files for ID: "x".'
        self.assertEqual(str(ctx.exception), msg)
        remove(join(self.tmpdir, 'x.faa'))
        remove(join(self.tmpdir, 'x.tsv'))

        # real files
        dir_ = join(self.datdir, 'align', 'bowtie2')

        obs = id2file_from_dir(dir_)
        exp = {f'S0{i}': f'S0{i}.sam.xz' for i in range(1, 6)}
        self.assertDictEqual(obs, exp)

        obs = id2file_from_dir(dir_, ext='.xz')
        exp = {f'S0{i}.sam': f'S0{i}.sam.xz' for i in range(1, 6)}
        self.assertDictEqual(obs, exp)

        obs = id2file_from_dir(dir_, ids=['S01', 'S02', 'S03'])
        exp = {f'S0{i}': f'S0{i}.sam.xz' for i in range(1, 4)}
        self.assertDictEqual(obs, exp)

    def test_id2file_from_map(self):
        # mock alignment files
        ids = ['a', 'b', 'c']
        fps = [join(self.tmpdir, '{}.fq'.format(x)) for x in ids]
        for fp in fps:
            open(fp, 'a').close()

        # mapping file with absolute paths
        mapfile = join(self.tmpdir, 'map.tmp')
        with open(mapfile, 'w') as f:
            for id_, fp in zip(ids, fps):
                f.write('{}\t{}\n'.format(id_, fp))
        obs = id2file_from_map(mapfile)
        exp = [pair for pair in zip(ids, fps)]
        self.assertListEqual(obs, exp)

        # relative paths to same directory
        with open(mapfile, 'w') as f:
            f.write('# This is a mapping file!\n\n')
            for id_ in ids:
                f.write('{}\t{}.fq\n'.format(id_, id_))
        obs = id2file_from_map(mapfile)
        self.assertListEqual(obs, exp)

        # not a mapping file
        with open(mapfile, 'w') as f:
            f.write('This is NOT a mapping file!\n\n')
        obs = id2file_from_map(mapfile)
        self.assertIsNone(obs)

        # incorrect filepaths starting from 1st line
        with open(mapfile, 'w') as f:
            for id_ in ids:
                f.write('{}\t{}.pdf\n'.format(id_, id_))
        obs = id2file_from_map(mapfile)
        self.assertIsNone(obs)

        # incorrect filepath(s) in middle of file
        with open(mapfile, 'w') as f:
            f.write('a\ta.fq\n')
            f.write('b\tb.fq\n')
            f.write('d\td.fq\n')
        with self.assertRaises(ValueError) as ctx:
            id2file_from_map(mapfile)
        self.assertEqual(str(ctx.exception), (
            'Alignment file "d.fq" does not exist.'))

        # empty mapping file
        open(mapfile, 'w').close()
        obs = id2file_from_map(mapfile)
        self.assertIsNone(obs)
        for fp in fps:
            remove(fp)

        # real files
        dir_ = join(self.datdir, 'align', 'bowtie2')
        exp = []
        with open(mapfile, 'w') as f:
            for i in range(1, 6):
                id_ = f'S0{i}'
                fp = join(dir_, f'{id_}.sam.xz')
                f.write(f'{id_}\t{fp}\n')
                exp.append((id_, fp))
        obs = id2file_from_map(mapfile)
        self.assertListEqual(obs, exp)
        remove(mapfile)

    def test_read_map(self):
        # simplest map
        tsv = ('R1	A',
               'R2	B',
               'R3	C')
        obs = read_map(tsv)
        exp = (('R1', 'A'), ('R2', 'B'), ('R3', 'C'))
        for obs_, exp_ in zip(read_map(tsv), exp):
            self.assertTupleEqual(obs_, exp_)

        # specify separator
        obs = read_map(('R1,A',), sep=',')
        self.assertTupleEqual(next(obs), ('R1', 'A'))

        # single column
        obs = tuple(read_map(('Hello.',)))
        self.assertTupleEqual(obs, ())

        # more than two columns and line is skipped
        obs = tuple(read_map(('Here\tit\tis.',)))
        self.assertTupleEqual(obs, ())

        # more than two columns and they are omitted
        obs = next(read_map(('Here\tit\tis.',), multi=False))
        exp = ('Here', 'it')
        self.assertTupleEqual(obs, exp)

        # multiple columns
        tsv = ('R1	A	B	C',
               'R2	D	E',
               'R3	F')
        obs = read_map(tsv, multi=True)
        exp = (('R1', ('A', 'B', 'C')), ('R2', ('D', 'E')), ('R3', ('F',)))
        self.assertTupleEqual(tuple(obs), exp)

        # single column with count
        tsv = ('R1	A',
               'R2	B:2',
               'R3	C:3')
        obs = read_map(tsv, count=True)
        exp = (('R1', ('A', 1)), ('R2', ('B', 2)), ('R3', ('C', 3)))
        self.assertTupleEqual(tuple(obs), exp)

        # multiple columns with count
        tsv = ('R1	A:3	B:2	C',
               'R2	D:4	E',
               'R3	F')
        obs = read_map(tsv, multi=True, count=True)
        exp = (('R1', (('A', 3), ('B', 2), ('C', 1))),
               ('R2', (('D', 4), ('E', 1))),
               ('R3', (('F', 1),)))
        self.assertTupleEqual(tuple(obs), exp)

    def test_write_readmap(self):
        # typical read map
        rmap = {'R1': 'G1',
                'R2': 'G2',
                'R3': {'G1': 1, 'G2': 2},
                'R4': {'G3': 3}}
        fp = join(self.tmpdir, 'readmap.tsv')
        with open(fp, 'w') as f:
            write_readmap(f, rmap)
        with open(fp, 'r') as f:
            obs = sorted(f.read().splitlines())
        exp = ['R1\tG1', 'R2\tG2', 'R3\tG2:2\tG1:1', 'R4\tG3:3']
        self.assertListEqual(obs, exp)

        # with name dict
        namedic = {'G1': 'Ecoli', 'G2': 'Strep', 'G3': 'Kleb'}
        with open(fp, 'w') as f:
            write_readmap(f, rmap, namedic)
        with open(fp, 'r') as f:
            obs = sorted(f.read().splitlines())
        exp = ['R1\tEcoli', 'R2\tStrep', 'R3\tStrep:2\tEcoli:1', 'R4\tKleb:3']
        self.assertListEqual(obs, exp)
        remove(fp)

    def test_write_table(self):
        # default mode
        data = {'S1': {'G1': 4, 'G2': 5, 'G3': 8},
                'S2': {'G1': 2, 'G4': 3, 'G5': 7},
                'S3': {'G2': 3, 'G5': 5}}
        fp = join(self.tmpdir, 'profile.tsv')
        with open(fp, 'w') as f:
            write_table(f, data)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3',
               'G1\t4\t2\t0',
               'G2\t5\t0\t3',
               'G3\t8\t0\t0',
               'G4\t0\t3\t0',
               'G5\t0\t7\t5']
        self.assertListEqual(obs, exp)

        # with sample Ids
        samples = ['S3', 'S1']
        with open(fp, 'w') as f:
            write_table(f, data, samples=samples)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS3\tS1',
               'G1\t0\t4',
               'G2\t3\t5',
               'G3\t0\t8',
               'G4\t0\t0',
               'G5\t5\t0']
        self.assertListEqual(obs, exp)

        # with taxon names
        namedic = {'G1': 'Actinobacteria',
                   'G2': 'Firmicutes',
                   'G3': 'Bacteroidetes',
                   'G4': 'Cyanobacteria'}
        with open(fp, 'w') as f:
            write_table(f, data, namedic=namedic)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3\tName',
               'G1\t4\t2\t0\tActinobacteria',
               'G2\t5\t0\t3\tFirmicutes',
               'G3\t8\t0\t0\tBacteroidetes',
               'G4\t0\t3\t0\tCyanobacteria',
               'G5\t0\t7\t5\t']
        self.assertListEqual(obs, exp)

        # with taxon names to replace Ids
        with open(fp, 'w') as f:
            write_table(f, data, namedic=namedic, name_as_id=True)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3',
               'Actinobacteria\t4\t2\t0',
               'Firmicutes\t5\t0\t3',
               'Bacteroidetes\t8\t0\t0',
               'Cyanobacteria\t0\t3\t0',
               'G5\t0\t7\t5']
        self.assertListEqual(obs, exp)

        # with ranks
        rankdic = {'G1': 'class', 'G2': 'phylum', 'G4': 'phylum'}
        with open(fp, 'w') as f:
            write_table(f, data, rankdic=rankdic)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3\tRank',
               'G1\t4\t2\t0\tclass',
               'G2\t5\t0\t3\tphylum',
               'G3\t8\t0\t0\t',
               'G4\t0\t3\t0\tphylum',
               'G5\t0\t7\t5\t']
        self.assertListEqual(obs, exp)

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
        with open(fp, 'w') as f:
            write_table(f, data, tree=tree)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3\tLineage',
               'G1\t4\t2\t0\t2;72;74',
               'G2\t5\t0\t3\t2;72',
               'G3\t8\t0\t0\t2;70',
               'G4\t0\t3\t0\t2;72',
               'G5\t0\t7\t5\t']
        self.assertListEqual(obs, exp)

        # with lineages and names as Ids
        namedic.update({
            '74': 'Actino', '72': 'Terra', '70': 'FCB', '2': 'Bacteria'})
        with open(fp, 'w') as f:
            write_table(f, data, tree=tree, namedic=namedic, name_as_id=True)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3\tLineage',
               'Actinobacteria\t4\t2\t0\tBacteria;Terra;Actino',
               'Firmicutes\t5\t0\t3\tBacteria;Terra',
               'Bacteroidetes\t8\t0\t0\tBacteria;FCB',
               'Cyanobacteria\t0\t3\t0\tBacteria;Terra',
               'G5\t0\t7\t5\t']
        self.assertListEqual(obs, exp)

        # with stratification
        sdata = {'S1': {('A', 'G1'): 4,
                        ('A', 'G2'): 5,
                        ('B', 'G1'): 8},
                 'S2': {('A', 'G1'): 2,
                        ('B', 'G1'): 3,
                        ('B', 'G2'): 7},
                 'S3': {('B', 'G3'): 3,
                        ('C', 'G2'): 5}}
        with open(fp, 'w') as f:
            write_table(f, sdata)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2\tS3',
               'A|G1\t4\t2\t0',
               'A|G2\t5\t0\t0',
               'B|G1\t8\t3\t0',
               'B|G2\t0\t7\t0',
               'B|G3\t0\t0\t3',
               'C|G2\t0\t0\t5']
        self.assertListEqual(obs, exp)

        # clean up
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
