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

from woltka.cli import gotu


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
            obs = [x.split('\t') for x in f.read().splitlines()]
        self.assertEqual(len(obs), 50)
        self.assertListEqual(obs[0], [
            '#FeatureID', 'S01', 'S02', 'S03', 'S04', 'S05'])
        self.assertListEqual(obs[1], [
            'G000006845', '0', '0', '10', '0', '0'])
        self.assertListEqual(obs[2], [
            'G000006925', '0', '6', '42', '0', '2'])
        remove(output_fp)


if __name__ == '__main__':
    main()
