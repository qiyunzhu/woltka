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

from click.testing import CliRunner

from woltka.cli import gotu, classify


class CliTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')
        self.runner = CliRunner()

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_gotu(self):
        output_fp = join(self.tmpdir, 'output.tsv')
        params = ['--input',  join(self.datdir, 'align', 'bowtie2'),
                  '--output', output_fp]
        res = self.runner.invoke(gotu, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1], 'Task completed.')
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'bowtie2.gotu.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)
        obs = [x.split('\t') for x in obs]
        self.assertEqual(len(obs), 50)
        self.assertListEqual(obs[0], [
            '#FeatureID', 'S01', 'S02', 'S03', 'S04', 'S05'])
        self.assertListEqual(obs[1], [
            'G000006845', '0', '0', '10', '0', '0'])
        self.assertListEqual(obs[2], [
            'G000006925', '0', '6', '42', '0', '2'])
        remove(output_fp)

    def test_classify(self):
        output_fp = join(self.tmpdir, 'output.tsv')

        # burst, classification at genus level
        params = ['--input',  join(self.datdir, 'align', 'burst'),
                  '--output', output_fp,
                  '--names',  join(self.datdir, 'taxonomy', 'names.dmp'),
                  '--nodes',  join(self.datdir, 'taxonomy', 'nodes.dmp'),
                  '--map',    join(self.datdir, 'taxonomy', 'g2tid.txt'),
                  '--rank',   'genus',
                  '--name-as-id']
        res = self.runner.invoke(classify, params)
        self.assertEqual(res.exit_code, 0)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'burst.genus.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)

        # bowtie2, free-rank classification
        params = ['--input',  join(self.datdir, 'align', 'bowtie2'),
                  '--output', output_fp,
                  '--nodes',  join(self.datdir, 'taxonomy', 'nodes.dmp'),
                  '--map',    join(self.datdir, 'taxonomy', 'g2tid.txt'),
                  '--rank',   'free',
                  '--no-subok']
        res = self.runner.invoke(classify, params)
        self.assertEqual(res.exit_code, 0)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'bowtie2.free.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)

        # bt2sho phylogeny-based classification
        params = ['--input',  join(self.datdir, 'align', 'bt2sho'),
                  '--output', output_fp,
                  '--newick', join(self.datdir, 'tree.nwk'),
                  '--rank',   'free']
        res = self.runner.invoke(classify, params)
        self.assertEqual(res.exit_code, 0)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'bt2sho.phylo.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)

        # burst, classification by GO process
        params = ['--input',  join(self.datdir, 'align', 'burst'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--coords', join(self.datdir, 'function', 'coords.txt.xz'),
                  '--map',    join(self.datdir, 'function', 'uniref.map.xz'),
                  '--map',    join(
                      self.datdir, 'function', 'go', 'process.tsv.xz'),
                  '--map-as-rank']
        res = self.runner.invoke(classify, params)
        self.assertEqual(res.exit_code, 0)
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'burst.process.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)
        remove(output_fp)


if __name__ == '__main__':
    main()
