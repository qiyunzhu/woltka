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

from woltka.cli import (
    cli, classify_cmd, normalize_cmd, filter_cmd, merge_cmd, collapse_cmd,
    coverage_cmd)


class CliTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')
        self.alndir = join(self.datdir, 'align')
        self.taxdir = join(self.datdir, 'taxonomy')
        self.fundir = join(self.datdir, 'function')
        self.outdir = join(self.datdir, 'output')
        self.runner = CliRunner()

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_cli(self):
        self.assertRaises(SystemExit, cli)

    def test_classify(self):
        output_fp = join(self.tmpdir, 'output.tsv')

        def _test_params(params, exp):
            res = self.runner.invoke(classify_cmd, params + ['--no-exe'])
            self.assertEqual(res.exit_code, 0)
            self.assertTrue(cmp(output_fp, join(self.outdir, exp)))

        # bowtie2, ogu
        params = ['--input',  join(self.alndir, 'bowtie2'),
                  '--output', output_fp]
        _test_params(params, 'bowtie2.ogu.tsv')

        # bowtie2, paired-end, free-rank classification
        params = ['--input',  join(self.alndir, 'bowtie2'),
                  '--output', output_fp,
                  '--nodes',  join(self.taxdir, 'nodes.dmp'),
                  '--map',    join(self.taxdir, 'taxid.map'),
                  '--rank',   'free']
        _test_params(params, 'bowtie2.free.tsv')

        # blastn, multiplexed, lineage-based, species-level classification
        params = ['--input',   join(self.alndir, 'blastn', 'mux.b6o.xz'),
                  '--output',  output_fp,
                  '--lineage', join(self.taxdir, 'lineages.txt'),
                  '--rank',    'species']
        _test_params(params, 'blastn.species.tsv')

        # burst, classification at genus level
        params = ['--input',  join(self.alndir, 'burst'),
                  '--output', output_fp,
                  '--outmap', self.tmpdir,
                  '--names',  join(self.taxdir, 'names.dmp'),
                  '--nodes',  join(self.taxdir, 'nodes.dmp'),
                  '--map',    join(self.taxdir, 'taxid.map'),
                  '--rank',   'genus',
                  '--name-as-id']
        _test_params(params, 'burst.genus.tsv')

        # check read maps
        for i in range(1, 6):
            outmap_fp = join(self.tmpdir, f'S0{i}.txt.gz')
            with gzip.open(outmap_fp, 'r') as f:
                obs = f.read()
            with gzip.open(join(self.outdir, 'burst.genus.map',
                                f'S0{i}.txt.gz'), 'r') as f:
                exp = f.read()
            self.assertEqual(obs, exp)
            remove(outmap_fp)

        # blastn, family-level classification, report relative abundance (%)
        params = ['--input',  join(self.alndir, 'blastn', 'mux.b6o.xz'),
                  '--output', output_fp,
                  '--names',  join(self.taxdir, 'names.dmp'),
                  '--nodes',  join(self.taxdir, 'nodes.dmp'),
                  '--map',    join(self.taxdir, 'taxid.map'),
                  '--rank',   'family',
                  '--name-as-id',
                  '--frac',
                  '--scale',  '100',
                  '--digits', 2]
        _test_params(params, 'blastn.family.percent.tsv')

        # bt2sho, order-level, genome size-normalized classification
        # unit: reads per million bases (of genome)
        params = ['--input',  join(self.alndir, 'bt2sho'),
                  '--output', output_fp,
                  '--names',  join(self.taxdir, 'names.dmp'),
                  '--nodes',  join(self.taxdir, 'nodes.dmp'),
                  '--map',    join(self.taxdir, 'taxid.map'),
                  '--rank',   'order',
                  '--sizes',  join(self.taxdir, 'length.map'),
                  '--scale',  '1M',
                  '--digits', 3]
        _test_params(params, 'bt2sho.order.cpm.tsv')

        # bt2sho phylogeny-based classification
        params = ['--input',  join(self.alndir, 'bt2sho'),
                  '--output', output_fp,
                  '--newick', join(self.datdir, 'tree.nwk'),
                  '--rank',   'free',
                  '--subok']
        _test_params(params, 'bt2sho.phylo.tsv')

        # burst, classification by GO process
        params = ['--input',  join(self.alndir, 'burst'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--coords', join(self.fundir, 'coords.txt.xz'),
                  '--map',    join(self.fundir, 'uniref', 'uniref.map.xz'),
                  '--map',    join(self.fundir, 'go', 'process.tsv.xz')]
        _test_params(params, 'burst.process.tsv')

        # burst, stratified genus/process classification
        params = ['--input',  join(self.alndir, 'burst'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--coords', join(self.fundir, 'coords.txt.xz'),
                  '--map',    join(self.fundir, 'uniref', 'uniref.map.xz'),
                  '--map',    join(self.fundir, 'go', 'process.tsv.xz'),
                  '--stratify', join(self.outdir, 'burst.genus.map')]
        _test_params(params, 'burst.genus.process.tsv')

        # bt2sho, classification by GO component, normalization by gene length
        # unit: reads per kilobase (RPK)
        params = ['--input',  join(self.alndir, 'bt2sho'),
                  '--output', output_fp,
                  '--rank',   'component',
                  '--coords', join(self.fundir, 'coords.txt.xz'),
                  '--map',    join(self.fundir, 'uniref', 'uniref.map.xz'),
                  '--map',    join(self.fundir, 'go', 'component.tsv.xz'),
                  '--sizes',  '.',
                  '--scale',  '1k',
                  '--digits', 3]
        _test_params(params, 'bt2sho.component.rpk.tsv')

        # simple map (from burst) against genes, classification at genus level
        params = ['--input',    join(self.alndir, 'burst', 'split'),
                  '--trim-sub', '_',
                  '--output',   output_fp,
                  '--rank',     'genus',
                  '--map',      join(self.taxdir, 'nucl', 'nucl2tid.txt'),
                  '--names',    join(self.taxdir, 'names.dmp'),
                  '--nodes',    join(self.taxdir, 'nodes.dmp'),
                  '--name-as-id']
        _test_params(params, 'split.genus.tsv')

        # simple map against genes, classification by GO process
        params = ['--input',  join(self.alndir, 'burst', 'split'),
                  '--output', output_fp,
                  '--rank',   'process',
                  '--map',    join(self.fundir, 'nucl', 'uniref.map.xz'),
                  '--map',    join(self.fundir, 'go', 'process.tsv.xz')]
        _test_params(params, 'split.process.tsv')

        remove(output_fp)

    def test_classify_cmd(self):
        input_fp = join(self.datdir, 'output', 'bowtie2.ogu.tsv')
        output_fp = join(self.tmpdir, 'tmp.tsv')
        sizes_fp = join(self.datdir, 'taxonomy', 'length.map')
        params = ['--input',  input_fp,
                  '--output', output_fp,
                  '--sizes',  sizes_fp,
                  '--scale',  '1M']
        res = self.runner.invoke(normalize_cmd, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1],
                         'Normalized profile written.')
        with open(output_fp, 'r') as f:
            obs = f.read().splitlines()
        self.assertIn('G000007325\t0\t12\t0\t0\t104', obs)
        self.assertIn('G000215745\t0\t0\t0\t358\t0', obs)
        remove(output_fp)

    def test_filter_cmd(self):
        input_fp = join(self.outdir, 'bowtie2.free.tsv')
        output_fp = join(self.tmpdir, 'output.tsv')
        params = ['--input',  input_fp,
                  '--output', output_fp,
                  '--min-percent', 1]
        res = self.runner.invoke(filter_cmd, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1],
                         'Filtered profile written.')
        with open(output_fp, 'r') as f:
            self.assertEqual(len(f.read().splitlines()), 30)
        remove(output_fp)

    def test_merge_cmd(self):
        input_fp1 = join(self.outdir, 'burst.process.tsv')
        input_fp2 = join(self.outdir, 'split.process.tsv')
        output_fp = join(self.tmpdir, 'output.tsv')
        params = ['--input',  input_fp1,
                  '--input',  input_fp2,
                  '--output', output_fp]
        res = self.runner.invoke(merge_cmd, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1],
                         'Merged profile written.')
        with open(output_fp, 'r') as f:
            self.assertEqual(len(f.read().splitlines()), 205)
        remove(output_fp)

    def test_collapse_cmd(self):
        input_fp = join(self.outdir, 'truth.uniref.tsv')
        map_fp = join(self.fundir, 'go', 'goslim.tsv.xz')
        output_fp = join(self.tmpdir, 'output.tsv')
        names_fp = join(self.fundir, 'go', 'name.txt.xz')
        params = ['--input',  input_fp,
                  '--map',    map_fp,
                  '--output', output_fp,
                  '--names',  names_fp,
                  '--divide']
        res = self.runner.invoke(collapse_cmd, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1],
                         'Collapsed profile written.')
        with open(output_fp, 'r') as f:
            self.assertEqual(len(f.read().splitlines()), 74)
        remove(output_fp)

    def test_coverage_cmd(self):
        input_fp = join(self.outdir, 'truth.metacyc.tsv')
        map_fp = join(self.fundir, 'metacyc', 'pathway_mbrs.txt')
        output_fp = join(self.tmpdir, 'output.tsv')
        names_fp = join(self.fundir, 'metacyc', 'pathway_name.txt')
        params = ['--input',  input_fp,
                  '--map',    map_fp,
                  '--output', output_fp,
                  '--names',  names_fp]
        res = self.runner.invoke(coverage_cmd, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output.splitlines()[-1],
                         'Coverage table written.')
        with open(output_fp, 'r') as f:
            self.assertEqual(len(f.read().splitlines()), 316)
        remove(output_fp)


if __name__ == '__main__':
    main()
