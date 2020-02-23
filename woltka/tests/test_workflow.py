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

from woltka.workflow import (
    parse_samples, build_align_proc, classify, workflow)


class WorkflowTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_parse_samples(self):
        # file (assuming demultiplexed)
        fp = join(self.tmpdir, 'input.fq')
        open(fp, 'a').close()
        obs = parse_samples(fp)
        self.assertIsNone(obs[0])
        self.assertListEqual(obs[1], [fp])
        self.assertTrue(obs[2])

        # confirm demux
        obs = parse_samples(fp, demux=True)
        self.assertIsNone(obs[0])
        self.assertListEqual(obs[1], [fp])
        self.assertTrue(obs[2])

        # force non-demux
        obs = parse_samples(fp, demux=False)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})
        self.assertFalse(obs[2])

        # provide correct extension
        obs = parse_samples(fp, ext='.fq', demux=False)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})

        # provide wrong extension
        with self.assertRaises(ValueError) as ctx:
            parse_samples(fp, ext='.fastq', demux=False)
        self.assertEqual(str(ctx.exception), (
            'Filepath and filename extension do not match.'))

        # provide Id list (no effect since demux)
        obs = parse_samples(fp, ids=['input'])
        self.assertListEqual(obs[0], ['input'])
        self.assertListEqual(obs[1], [fp])

        # provide correct Id list
        obs = parse_samples(fp, ids=['input'], demux=False)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})

        # provide wrong Id list
        with self.assertRaises(ValueError) as ctx:
            parse_samples(fp, ids=['hello'], demux=False)
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))

        # directory (assuming per-sample)
        obs = parse_samples(self.tmpdir)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})
        self.assertFalse(obs[2])

        # directory empty
        remove(fp)
        with self.assertRaises(ValueError) as ctx:
            parse_samples(self.tmpdir)
        self.assertEqual(str(ctx.exception), (
            'No valid file found in directory.'))

        # multiple files
        for i in range(1, 4):
            open(join(self.tmpdir, f'S{i}.sam'), 'a').close()
        obs = parse_samples(self.tmpdir)[1]
        exp = {join(self.tmpdir, f'S{i}.sam'): f'S{i}' for i in range(1, 4)}
        self.assertDictEqual(obs, exp)

        # add an irrelevant file
        fp = join(self.tmpdir, 'readme.txt')
        open(fp, 'a').close()
        obs = parse_samples(self.tmpdir)[1]
        exp[fp] = 'readme'
        self.assertDictEqual(obs, exp)

        # specify extension to target alignment files only
        del exp[fp]
        obs = parse_samples(self.tmpdir, ext='.sam')[1]
        self.assertDictEqual(obs, exp)

        # specify sample Ids
        obs = parse_samples(self.tmpdir, ids=['S1', 'S2', 'S3'])[1]
        self.assertDictEqual(obs, exp)

        # some samples are not found
        remove(fp)
        with self.assertRaises(ValueError) as ctx:
            parse_samples(self.tmpdir, ids=['S1', 'S2', 'S4'])
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))

        # force demux
        obs = parse_samples(self.tmpdir, demux=True)
        self.assertIsNone(obs[0])
        self.assertTrue(obs[2])
        exp = [join(self.tmpdir, f'S{i}.sam') for i in range(1, 4)]
        self.assertListEqual(obs[1], exp)

        # sample Ids are ignored when demux
        obs = parse_samples(self.tmpdir, ids=['S1', 'S2', 'S4'], demux=True)
        self.assertListEqual(obs[0], ['S1', 'S2', 'S4'])
        exp = [join(self.tmpdir, f'S{i}.sam') for i in range(1, 4)]
        self.assertListEqual(obs[1], exp)

        # not a valid path
        with self.assertRaises(ValueError) as ctx:
            parse_samples('im/not/path')
        self.assertEqual(str(ctx.exception), (
            '"im/not/path" is not a valid file or directory.'))

    def test_workflow(self):
        # simplest gotu workflow
        input_path = join(self.datdir, 'align', 'bowtie2')
        output_path = join(self.tmpdir, 'tmp.tsv')
        obs = workflow(input_path, output_path)['none']
        self.assertEqual(obs['S01']['G000011545'], 48)
        self.assertNotIn('G000007145', obs['S02'])
        self.assertEqual(obs['S03']['G000009345'], 640)
        with open(output_path, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'bowtie2.gotu.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)
        remove(output_path)

    def test_classify(self):
        # simplest gotu workflow
        input_path = join(self.datdir, 'align', 'bowtie2')
        samples, files, demux = parse_samples(input_path)
        proc = build_align_proc()
        ranks = ['none']
        obs = classify(proc, files, samples=samples, demux=demux,
                       ranks=ranks)['none']
        self.assertEqual(obs['S01']['G000011545'], 48)
        self.assertNotIn('G000007145', obs['S02'])
        self.assertEqual(obs['S03']['G000009345'], 640)
        self.assertEqual(obs['S04']['G000240185'], 4)
        self.assertEqual(obs['S05']['G000191145'], 10)


if __name__ == '__main__':
    main()
