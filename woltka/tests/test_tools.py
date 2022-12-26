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
from shutil import rmtree, copy
from tempfile import mkdtemp

from woltka.tools import (
    normalize_wf, filter_wf, merge_wf, collapse_wf, coverage_wf)


class ToolsTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_normalize_wf(self):
        input_fp = join(self.datdir, 'output', 'bowtie2.ogu.tsv')
        output_fp = join(self.tmpdir, 'tmp.tsv')

        # relative abundance
        normalize_wf(input_fp, output_fp, digits=3)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertIn('G000007325\t0.0\t0.013\t0.0\t0.0\t0.113', obs)
        self.assertIn('G000215745\t0.0\t0.0\t0.0\t0.944\t0.0', obs)

        # by genome size (cpm)
        sizes_fp = join(self.datdir, 'taxonomy', 'length.map')
        normalize_wf(input_fp, output_fp, sizes_fp, scale='1M')
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertIn('G000007325\t0\t12\t0\t0\t104', obs)
        self.assertIn('G000215745\t0\t0\t0\t358\t0', obs)

        # bad scale param
        with self.assertRaises(SystemExit) as ctx:
            normalize_wf(input_fp, output_fp, scale='Hi!')
        errmsg = '"Hi!" is not a valid scale factor.'
        self.assertEqual(str(ctx.exception), errmsg)

        # by gene size (from coords)
        input_fp = join(self.datdir, 'output', 'bowtie2.orf.tsv')
        sizes_fp = join(self.datdir, 'function', 'coords.txt.xz')
        normalize_wf(input_fp, output_fp, sizes_fp, scale='1k', digits=3)
        with open(output_fp, 'r') as f:
            for line in f:
                gene, _, datum = line.rstrip().partition('\t')
                if gene == 'G000006925_52':
                    self.assertEqual(datum, '0.0\t0.0\t0.849\t0.0\t0.0')
                if gene == 'G000008165_4079':
                    self.assertEqual(datum, '0.0\t0.0\t1.178\t0.0\t0.0')
                if gene == 'G000014525_1009':
                    self.assertEqual(datum, '0.0\t0.0\t0.0\t0.0\t0.887')

        # missing sizes
        sizes_fp = join(self.datdir, 'tree.nwk')
        with self.assertRaises(SystemExit) as ctx:
            normalize_wf(input_fp, output_fp, sizes_fp)
        errmsg = 'One or more features are not found in the size map.'
        self.assertEqual(str(ctx.exception), errmsg)

        # empty size file
        open(output_fp, 'w').close()
        with self.assertRaises(SystemExit) as ctx:
            normalize_wf(input_fp, output_fp, output_fp)
        errmsg = 'Size map file is empty or unreadable.'
        self.assertEqual(str(ctx.exception), errmsg)
        remove(output_fp)

    def test_filter_wf(self):
        input_fp = join(self.datdir, 'output', 'blastn.species.tsv')
        output_fp = join(self.tmpdir, 'tmp.tsv')

        # filter by count
        filter_wf(input_fp, output_fp, min_count=5)
        with open(output_fp, 'r') as f:
            n = len(f.read().splitlines())
        self.assertEqual(n, 44)

        # filter by percentage
        filter_wf(input_fp, output_fp, min_percent=1)
        with open(output_fp, 'r') as f:
            n = len(f.read().splitlines())
        self.assertEqual(n, 31)
        remove(output_fp)

        # no threshold specified
        with self.assertRaises(SystemExit) as ctx:
            filter_wf(input_fp, output_fp)
        errmsg = ('Please specify either minimum count or minimum percentage '
                  'threshold.')
        self.assertEqual(str(ctx.exception), errmsg)

        # zero is not a valid threshold
        with self.assertRaises(SystemExit) as ctx:
            filter_wf(input_fp, output_fp, min_count=0, min_percent=0)
        self.assertEqual(str(ctx.exception), errmsg)

        # both thresholds specified
        with self.assertRaises(SystemExit) as ctx:
            filter_wf(input_fp, output_fp, min_count=10, min_percent=10)
        errmsg = ('Only one of minimum count or minimum percentage thresholds '
                  'can be specified.')
        self.assertEqual(str(ctx.exception), errmsg)

        # percentage threshold too large
        with self.assertRaises(SystemExit) as ctx:
            filter_wf(input_fp, output_fp, min_percent=120)
        errmsg = 'Minimum percentage threshold must be below 100.'
        self.assertEqual(str(ctx.exception), errmsg)

    def test_merge_wf(self):
        # files
        input_fp1 = join(self.datdir, 'output', 'burst.process.tsv')
        input_fp2 = join(self.datdir, 'output', 'split.process.tsv')
        output_fp = join(self.tmpdir, 'tmp.tsv')
        merge_wf([input_fp1, input_fp2], output_fp)
        with open(output_fp, 'r') as f:
            n = len(f.read().splitlines())
        self.assertEqual(n, 205)
        remove(output_fp)

        # directory
        input_fp = join(self.tmpdir, 'indir')
        makedirs(input_fp)

        # only one profile
        copy(input_fp1, join(input_fp, 'test1.tsv'))
        with self.assertRaises(SystemExit) as ctx:
            merge_wf([input_fp], output_fp)
        errmsg = 'Please provide two or more profiles.'
        self.assertEqual(str(ctx.exception), errmsg)

        # two profiles
        copy(input_fp2, join(input_fp, 'test2.tsv'))
        merge_wf([input_fp], output_fp)
        with open(output_fp, 'r') as f:
            n = len(f.read().splitlines())
        self.assertEqual(n, 205)
        remove(output_fp)

        # invalid profile
        copy(join(self.datdir, 'tree.nwk'), join(input_fp, 'test3.tsv'))
        with self.assertRaises(SystemExit) as ctx:
            merge_wf([input_fp], output_fp)
        errmsg = 'Cannot parse test3.tsv as a profile.'
        self.assertEqual(str(ctx.exception), errmsg)

        # keep original decimal precision
        input_fp1 = join(self.tmpdir, '1st.tsv')
        with open(input_fp1, 'w') as f:
            f.write('#FeatureID\tS01\nG01\t1.1')
        input_fp2 = join(self.tmpdir, '2nd.tsv')
        with open(input_fp2, 'w') as f:
            f.write('#FeatureID\tS01\nG01\t2.2')
        output_fp = join(self.tmpdir, 'tmp.tsv')
        merge_wf([input_fp1, input_fp2], output_fp)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS01', 'G01\t3.3']
        self.assertListEqual(obs, exp)
        remove(output_fp)
        remove(input_fp1)
        remove(input_fp2)

    def test_collapse_wf(self):
        input_fp = join(self.datdir, 'output', 'truth.gene.tsv')
        output_fp = join(self.tmpdir, 'tmp.tsv')
        map_fp = join(self.datdir, 'function', 'nucl', 'uniref.map.xz')

        # one-to-one
        collapse_wf(input_fp, output_fp, map_fp)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertEqual(len(obs), 2489)
        for line in obs:
            if line.startswith('A0A069B1G8'):
                self.assertEqual(line[11:], '1\t0\t0\t0\t0')
            if line.startswith('P0AGG9\t'):
                self.assertEqual(line[7:], '0\t0\t0\t2\t0')

        # one-to-many
        input_fp = join(self.datdir, 'output', 'truth.uniref.tsv')
        map_fp = join(self.datdir, 'function', 'go', 'goslim.tsv.xz')
        names_fp = join(self.datdir, 'function', 'go', 'name.txt.xz')
        collapse_wf(input_fp, output_fp, map_fp, names_fp=names_fp)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertEqual(len(obs), 84)
        for line in obs:
            if line.startswith('GO:0005737'):
                self.assertEqual(line[11:], '59\t55\t62\t22\t47\tcytoplasm')
            if line.startswith('GO:0004518'):
                self.assertEqual(line[11:], '3\t7\t9\t1\t4\tnuclease activity')
        remove(output_fp)

        # wrong mapping file
        map_fp = join(self.datdir, 'tree.nwk')
        with self.assertRaises(SystemExit) as ctx:
            collapse_wf(input_fp, output_fp, map_fp, divide=True)
        errmsg = 'No source-target mapping is found in tree.nwk.'
        self.assertEqual(str(ctx.exception), errmsg)

        # stratified profile
        input_fp = join(self.datdir, 'output', 'burst.genus.process.tsv')
        map_fp = join(self.datdir, 'function', 'go', 'go2slim.map.xz')
        collapse_wf(input_fp, output_fp, map_fp, field=2)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertEqual(len(obs), 279)
        for line in obs:
            if line.startswith('Klebsiella|GO:0008150'):
                self.assertEqual(line[22:], '2\t47\t0\t87\t0')
            if line.startswith('Streptococcus|GO:0008150'):
                self.assertEqual(line[25:], '0\t2\t9\t3\t0')
        remove(output_fp)

        # nested profile
        input_fp = join(self.datdir, 'output', 'bowtie2.orf.tsv')
        collapse_wf(input_fp, output_fp, field=1, sep='_', nested=True)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertEqual(len(obs), 50)
        for line in obs:
            if line.startswith('G000011705'):
                self.assertEqual(line[11:], '104\t0\t29\t0\t0')
            if line.startswith('G000240185'):
                self.assertEqual(line[11:], '2\t608\t2\t1\t0')

    def test_coverage_wf(self):
        input_fp = join(self.datdir, 'output', 'truth.metacyc.tsv')
        map_fp = join(self.datdir, 'function', 'metacyc', 'pathway_mbrs.txt')
        output_fp = join(self.tmpdir, 'tmp.tsv')

        # default behavior
        coverage_wf(input_fp, map_fp, output_fp)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertEqual(len(obs), 316)
        for line in obs:
            if line.startswith('ARABCAT-PWY\t'):
                self.assertListEqual(line.split('\t')[1:], [
                    '20.0', '80.0', '100.0', '80.0', '100.0'])

        # parameters
        names_fp = join(self.datdir, 'function', 'metacyc', 'pathway_name.txt')
        coverage_wf(input_fp, map_fp, output_fp, threshold=80,
                    names_fp=names_fp)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertEqual(len(obs), 86)
        for line in obs:
            if line.startswith('ARABCAT-PWY\t'):
                self.assertListEqual(line.split('\t')[1:], [
                    '0', '1', '1', '1', '1', 'L-arabinose degradation I'])
        remove(output_fp)

        # wrong mapping file
        map_fp = join(self.datdir, 'tree.nwk')
        with self.assertRaises(SystemExit) as ctx:
            coverage_wf(input_fp, map_fp, output_fp, count=True)
        errmsg = 'No group membership is found in tree.nwk.'
        self.assertEqual(str(ctx.exception), errmsg)


if __name__ == '__main__':
    main()
