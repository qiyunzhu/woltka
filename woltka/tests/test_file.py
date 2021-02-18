#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove, makedirs
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp
import gzip
import bz2

from woltka.file import (
    openzip, readzip, file2stem, path2stem, stem2rank, read_ids,
    id2file_from_dir, id2file_from_map, read_map_uniq, read_map_1st,
    read_map_all, read_map_many, write_readmap)


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

    def test_readzip(self):
        text = 'Hello World!'

        # read regular file
        fp = join(self.tmpdir, 'test.txt')
        with open(fp, 'w') as f:
            f.write(text)
        with readzip(fp) as f:
            obs = f.read()
        self.assertEqual(obs, text)

        # read gzip file using builti-in Python module
        fpz = join(self.tmpdir, 'test.txt.gz')
        with gzip.open(fpz, 'wb') as f:
            f.write(text.encode())
        with readzip(fpz) as f:
            obs = f.read()
        self.assertEqual(obs, text)

        # read gzip file using auto-identified gzip program
        zippers = {}
        with readzip(fpz, zippers) as f:
            obs = f.read()
        self.assertEqual(obs, text)
        self.assertDictEqual(zippers, {'gzip': True})

        # read gzip file using already-identified gzip program
        with readzip(fpz, zippers) as f:
            obs = f.read()
        self.assertEqual(obs, text)
        self.assertTrue(zippers['gzip'])

        # disable gzip program
        zippers['gzip'] = False
        with readzip(fpz, zippers) as f:
            obs = f.read()
        self.assertEqual(obs, text)
        self.assertFalse(zippers['gzip'])
        remove(fpz)

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

    def test_stem2rank(self):
        self.assertEqual(stem2rank('/path/to/db'), 'db')
        self.assertEqual(stem2rank('gene.map.gz'), 'gene')
        self.assertEqual(stem2rank('accession2taxid.txt.bz2'), 'taxid')
        self.assertEqual(stem2rank('compound_2_reaction.txt'), 'reaction')
        self.assertEqual(stem2rank('strain-to-species.map'), 'species')
        self.assertEqual(stem2rank('you_and_i.txt'), 'you_and_i')

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

        # skip subdirectory and only consider files
        sdir = join(self.tmpdir, 'im_dir')
        makedirs(sdir)
        obs = id2file_from_dir(self.tmpdir)
        self.assertDictEqual(obs, exp)
        rmtree(sdir)
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

    def test_read_map_uniq(self):
        # simplest map (two columns)
        tsv = ('R1	A',
               'R2	B',
               'R3	C')
        exp = (('R1', 'A'), ('R2', 'B'), ('R3', 'C'))
        for obs_, exp_ in zip(read_map_uniq(tsv), exp):
            self.assertTupleEqual(obs_, exp_)

        # specify separator
        obs = read_map_uniq(('R1,A',), sep=',')
        self.assertTupleEqual(next(obs), ('R1', 'A'))

        # single column
        obs = tuple(read_map_uniq(('Hello.',)))
        self.assertTupleEqual(obs, ())

        # more than two columns
        obs = tuple(read_map_uniq(('Here\tit\tis.',)))
        self.assertTupleEqual(obs, ())

    def test_read_map_1st(self):
        # simplest map (two columns)
        tsv = ('R1	A',
               'R2	B',
               'R3	C')
        exp = (('R1', 'A'), ('R2', 'B'), ('R3', 'C'))
        for obs_, exp_ in zip(read_map_1st(tsv), exp):
            self.assertTupleEqual(obs_, exp_)

        # specify separator
        obs = read_map_1st(('R1,A',), sep=',')
        self.assertTupleEqual(next(obs), ('R1', 'A'))

        # single column
        obs = tuple(read_map_1st(('Hello.',)))
        self.assertTupleEqual(obs, ())

        # more than two columns
        obs = read_map_1st(('Here\tit\tis.',))
        self.assertTupleEqual(next(obs), ('Here', 'it'))

    def test_read_map_all(self):
        # simple map (two or more columns)
        tsv = ('R1	A',
               'R2	A	B',
               'R3	C	A	B')
        exp = (('R1', ['A']), ('R2', ['A', 'B']), ('R3', ['C', 'A', 'B']))
        for obs_, exp_ in zip(read_map_all(tsv), exp):
            self.assertTupleEqual(obs_, exp_)

        # specify separator
        obs = read_map_all(('R1,A,B,C',), sep=',')
        self.assertTupleEqual(next(obs), ('R1', ['A', 'B', 'C']))

        # single column
        obs = tuple(read_map_all(('Hello.',)))
        self.assertTupleEqual(obs, ())

    def test_read_map_many(self):
        # one-to-one
        tsv = ('R1	A',
               'R2	B',
               'R3	C')
        obs = read_map_many(tsv)
        exp = {'R1': ['A'], 'R2': ['B'], 'R3': ['C']}
        self.assertDictEqual(obs, exp)

        # many-to-one, with separator
        tsv = ('R1,A',
               'R2,B',
               'R3,B',
               'R4,C')
        obs = read_map_many(tsv, sep=',')
        exp = {'R1': ['A'], 'R2': ['B'], 'R3': ['B'], 'R4': ['C']}
        self.assertDictEqual(obs, exp)

        # one-to-many
        tsv = ('R1	A',
               'R2	A	B',
               'R3	B	C	D',
               'R4	A	D')
        obs = read_map_many(tsv)
        exp = {'R1': ['A'], 'R2': ['A', 'B'], 'R3': ['B', 'C', 'D'],
               'R4': ['A', 'D']}
        self.assertDictEqual(obs, exp)

        # many-to-many
        tsv = ('R1	A',
               'R1	B',
               'R2	B',
               'R3	A',
               'R3	C',
               'R4	C')
        obs = read_map_many(tsv)
        exp = {'R1': ['A', 'B'], 'R2': ['B'], 'R3': ['A', 'C'], 'R4': ['C']}
        self.assertDictEqual(obs, exp)

        # no relationship
        tsv = ('Hello.')
        obs = read_map_many(tsv)
        self.assertDictEqual(obs, {})

    def test_write_readmap(self):
        # typical read map
        rmap = {'R1': 'G1',
                'R2': 'G2',
                'R3': {'G1': 1, 'G2': 2},
                'R4': {'G3': 3}}
        fp = join(self.tmpdir, 'readmap.tsv')
        with open(fp, 'w') as f:
            write_readmap(f, rmap.keys(), rmap.values())
        with open(fp, 'r') as f:
            obs = sorted(f.read().splitlines())
        exp = ['R1\tG1', 'R2\tG2', 'R3\tG2:2\tG1:1', 'R4\tG3:3']
        self.assertListEqual(obs, exp)

        # with name dict
        namedic = {'G1': 'Ecoli', 'G2': 'Strep', 'G3': 'Kleb'}
        with open(fp, 'w') as f:
            write_readmap(f, rmap.keys(), rmap.values(), namedic)
        with open(fp, 'r') as f:
            obs = sorted(f.read().splitlines())
        exp = ['R1\tEcoli', 'R2\tStrep', 'R3\tStrep:2\tEcoli:1', 'R4\tKleb:3']
        self.assertListEqual(obs, exp)
        remove(fp)


if __name__ == '__main__':
    main()
