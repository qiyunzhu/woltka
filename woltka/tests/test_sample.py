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

from woltka.sample import read_ids, match_sample_file, id2file_map


class SampleTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_read_ids(self):
        # simple list
        ids = ['a', 'b', 'c']
        obs = read_ids(ids)
        self.assertListEqual(obs, ids)

        # metadata table
        lines = ['#SampleID\tHeight\tWeight',
                 'S01\t120.5\t56.8',
                 'S02\t114.2\t52.5',
                 'S03\t117.8\t55.0']
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

    def test_match_sample_file(self):
        # file (assuming demultiplexed)
        fp = join(self.tmpdir, 'input.fq')
        open(fp, 'a').close()
        obs0, obs1 = match_sample_file(fp)
        self.assertTrue(obs0)
        self.assertListEqual(obs1, [fp])

        # confirm demux
        obs0, obs1 = match_sample_file(fp, demux=True)
        self.assertTrue(obs0)
        self.assertListEqual(obs1, [fp])

        # force non-demux
        obs0, obs1 = match_sample_file(fp, demux=False)
        self.assertFalse(obs0)
        self.assertDictEqual(obs1, {fp: 'input'})

        # provide correct extension
        obs = match_sample_file(fp, demux=False, ext='.fq')[1]
        self.assertDictEqual(obs1, {fp: 'input'})

        # provide wrong extension
        with self.assertRaises(ValueError) as ctx:
            match_sample_file(fp, demux=False, ext='.fastq')
        self.assertEqual(str(ctx.exception), (
            'Filepath and filename extension do not match.'))

        # provide Id list (no effect since demux)
        obs = match_sample_file(fp, samples=['input'])[1]
        self.assertListEqual(obs, [fp])

        # provide correct Id list
        obs = match_sample_file(fp, demux=False, samples=['input'])[1]
        self.assertDictEqual(obs, {fp: 'input'})

        # provide wrong Id list
        with self.assertRaises(ValueError) as ctx:
            match_sample_file(fp, demux=False, samples=['hello'])
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))

        # directory (assuming per-sample)
        obs0, obs1 = match_sample_file(self.tmpdir)
        self.assertFalse(obs0)
        self.assertDictEqual(obs1, {fp: 'input'})

        # directory empty
        remove(fp)
        with self.assertRaises(ValueError) as ctx:
            match_sample_file(self.tmpdir)
        self.assertEqual(str(ctx.exception), (
            'No valid file found in directory.'))

        # multiple files
        for i in range(1, 4):
            open(join(self.tmpdir, f'S{i}.sam'), 'a').close()
        obs = match_sample_file(self.tmpdir)[1]
        exp = {join(self.tmpdir, f'S{i}.sam'): f'S{i}' for i in range(1, 4)}
        self.assertDictEqual(obs, exp)

        # add an irrelevant file
        fp = join(self.tmpdir, 'readme.txt')
        open(fp, 'a').close()
        obs = match_sample_file(self.tmpdir)[1]
        exp[fp] = 'readme'
        self.assertDictEqual(obs, exp)

        # specify extension to target alignment files only
        del exp[fp]
        obs = match_sample_file(self.tmpdir, ext='.sam')[1]
        self.assertDictEqual(obs, exp)

        # specify sample Ids
        obs = match_sample_file(self.tmpdir, samples=['S1', 'S2', 'S3'])[1]
        self.assertDictEqual(obs, exp)

        # some samples are not found
        remove(fp)
        with self.assertRaises(ValueError) as ctx:
            match_sample_file(self.tmpdir, samples=['S1', 'S2', 'S4'])
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))

        # force demux
        obs0, obs1 = match_sample_file(self.tmpdir, demux=True)
        self.assertTrue(obs0)
        exp1 = [join(self.tmpdir, f'S{i}.sam') for i in range(1, 4)]
        self.assertListEqual(obs1, exp1)

        # sample Ids are ignored when demux
        obs = match_sample_file(self.tmpdir, demux=True, samples=[
            'S1', 'S2', 'S4'])[1]
        exp = [join(self.tmpdir, f'S{i}.sam') for i in range(1, 4)]
        self.assertListEqual(obs, exp)

        # not a valid path
        with self.assertRaises(ValueError) as ctx:
            match_sample_file('im/not/path')
        self.assertEqual(str(ctx.exception), (
            '"im/not/path" is not a valid file or directory.'))

    def test_id2file_map(self):
        ids = ['a', 'b', 'c']

        # regular Fasta files
        for id_ in ids:
            open(join(self.tmpdir, '{}.faa'.format(id_)), 'a').close()
        obs = id2file_map(self.tmpdir)
        exp = {x: '{}.faa'.format(x) for x in ids}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, '{}.faa'.format(id_)))

        # gzipped Fasta files
        for id_ in ids:
            open(join(self.tmpdir, '{}.faa.gz'.format(id_)), 'a').close()
        obs = id2file_map(self.tmpdir)
        exp = {x: '{}.faa.gz'.format(x) for x in ids}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, '{}.faa.gz'.format(id_)))

        # user-defined extension filename
        for id_ in ids:
            open(join(self.tmpdir, '{}.faa'.format(id_)), 'a').close()
        open(join(self.tmpdir, 'readme.txt'), 'a').close()
        open(join(self.tmpdir, 'taxdump'), 'a').close()
        obs = id2file_map(self.tmpdir, ext='.faa')
        exp = {x: '{}.faa'.format(x) for x in ids}
        self.assertDictEqual(obs, exp)
        remove(join(self.tmpdir, 'readme.txt'))
        remove(join(self.tmpdir, 'taxdump'))

        # user-defined Id list
        obs = id2file_map(self.tmpdir, ids=['a', 'b'])
        exp = {x: '{}.faa'.format(x) for x in ['a', 'b']}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, '{}.faa'.format(id_)))

        # duplicated Ids
        open(join(self.tmpdir, 'x.faa'), 'a').close()
        open(join(self.tmpdir, 'x.tsv'), 'a').close()
        with self.assertRaises(ValueError) as ctx:
            id2file_map(self.tmpdir)
        msg = 'Ambiguous files for ID: "x".'
        self.assertEqual(str(ctx.exception), msg)
        remove(join(self.tmpdir, 'x.faa'))
        remove(join(self.tmpdir, 'x.tsv'))

        # read files
        dir_ = join(self.datdir, 'align', 'bowtie2')

        obs = id2file_map(dir_)
        exp = {f'S0{i}': f'S0{i}.sam.xz' for i in range(1, 6)}
        self.assertDictEqual(obs, exp)

        obs = id2file_map(dir_, ext='.xz')
        exp = {f'S0{i}.sam': f'S0{i}.sam.xz' for i in range(1, 6)}
        self.assertDictEqual(obs, exp)

        obs = id2file_map(dir_, ids=['S01', 'S02', 'S03'])
        exp = {f'S0{i}': f'S0{i}.sam.xz' for i in range(1, 4)}
        self.assertDictEqual(obs, exp)


if __name__ == '__main__':
    main()
