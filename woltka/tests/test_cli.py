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
from filecmp import cmp
import gzip

from click.testing import CliRunner

from woltka.cli import cli, gotu_cmd, classify_cmd, filter_cmd


class CliTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')
        self.runner = CliRunner()

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_cli(self):
        self.assertRaises(SystemExit, cli)

    def test_gotu(self):
        output_fp = join(self.tmpdir, 'output.tsv')
        params = ['--input',  join(self.datdir, 'align', 'bowtie2'),
                  '--output', output_fp]
        res = self.runner.invoke(gotu_cmd, params)
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

        def _test_params(params, exp):
            res = self.runner.invoke(classify_cmd, params)
            self.assertEqual(res.exit_code, 0)
            self.assertTrue(cmp(output_fp, join(self.datdir, 'output', exp)))

        # bowtie2, paired-end, free-rank classification
        params = ['--input',  join(self.datdir, 'align', 'bowtie2'),
                  '--output', output_fp,
                  '--nodes',  join(self.datdir, 'taxonomy', 'nodes.dmp'),
                  '--map',    join(self.datdir, 'taxonomy', 'g2tid.txt'),
                  '--rank',   'free',
                  '--no-subok']
        _test_params(params, 'bowtie2.free.tsv')

        # blastn, multiplexed, lineage-based, species-level classification
        params = ['--input',   join(
                      self.datdir, 'align', 'blastn', 'mux.b6o.xz'),
                  '--output',  output_fp,
                  '--lineage', join(self.datdir, 'taxonomy', 'lineage.txt'),
                  '--rank',    'species']
        _test_params(params, 'blastn.species.tsv')

        # burst, classification at genus level
        params = ['--input',  join(self.datdir, 'align', 'burst'),
                  '--output', output_fp,
                  '--outmap', self.tmpdir,
                  '--names',  join(self.datdir, 'taxonomy', 'names.dmp'),
                  '--nodes',  join(self.datdir, 'taxonomy', 'nodes.dmp'),
                  '--map',    join(self.datdir, 'taxonomy', 'g2tid.txt'),
                  '--rank',   'genus',
                  '--name-as-id']
        _test_params(params, 'burst.genus.tsv')

        # check read maps
        for i in range(1, 6):
            outmap_fp = join(self.tmpdir, f'S0{i}.txt.gz')
            with gzip.open(outmap_fp, 'r') as f:
                obs = f.read()
            with gzip.open(join(self.datdir, 'output', 'burst.genus.map',
                                f'S0{i}.txt.gz'), 'r') as f:
                exp = f.read()
            self.assertEqual(obs, exp)
            remove(outmap_fp)

        # bt2sho phylogeny-based classification
        params = ['--input',  join(self.datdir, 'align', 'bt2sho'),
                  '--output', output_fp,
                  '--newick', join(self.datdir, 'tree.nwk'),
                  '--rank',   'free']
        _test_params(params, 'bt2sho.phylo.tsv')

        # burst, classification by GO process
        params = ['--input',  join(self.datdir, 'align', 'burst'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--coords', join(self.datdir, 'function', 'coords.txt.xz'),
                  '--map',    join(self.datdir, 'function', 'uniref.map.xz'),
                  '--map',    join(
                      self.datdir, 'function', 'go', 'process.tsv.xz'),
                  '--map-as-rank']
        _test_params(params, 'burst.process.tsv')

        # burst, stratified genus/process classification
        params = ['--input',  join(self.datdir, 'align', 'burst'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--coords', join(self.datdir, 'function', 'coords.txt.xz'),
                  '--map',    join(self.datdir, 'function', 'uniref.map.xz'),
                  '--map',    join(
                      self.datdir, 'function', 'go', 'process.tsv.xz'),
                  '--map-as-rank',
                  '--stratify', join(self.datdir, 'output', 'burst.genus.map')]
        _test_params(params, 'burst.genus.process.tsv')

        # simple map (from burst) against genes, classification at genus level
        params = ['--input',  join(self.datdir, 'align', 'burst', 'split'),
                  '--output', output_fp,
                  '--rank',   'genus',
                  '--map',    join(
                      self.datdir, 'taxonomy', 'nucl', 'nucl2tid.txt'),
                  '--names',  join(self.datdir, 'taxonomy', 'names.dmp'),
                  '--nodes',  join(self.datdir, 'taxonomy', 'nodes.dmp'),
                  '--name-as-id',
                  '--deidx']
        _test_params(params, 'split.genus.tsv')

        # simple map against genes, classification by GO process
        params = ['--input',  join(self.datdir, 'align', 'burst', 'split'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--map',    join(
                      self.datdir, 'function', 'nucl', 'uniref.map.xz'),
                  '--map',    join(
                      self.datdir, 'function', 'go', 'process.tsv.xz'),
                  '--map-as-rank']
        _test_params(params, 'split.process.tsv')

        remove(output_fp)

    def test_filter_table(self):
        input_fp = join(self.datdir, 'output', 'bowtie2.free.tsv')
        output_fp = join(self.tmpdir, 'output.tsv')
        params = ['--input', input_fp,
                  '--output', output_fp,
                  '--min-percent', 1]
        res = self.runner.invoke(filter_cmd, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1],
                         'Filtered profile written.')
        with open(output_fp, 'r') as f:
            self.assertEqual(len(f.read().splitlines()) - 1, 29)
        remove(output_fp)


if __name__ == '__main__':
    main()
