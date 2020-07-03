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

from woltka.tools import filter_wf


class ToolsTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_filter_wf(self):
        input_fp = join(self.datdir, 'output', 'blastn.species.tsv')
        output_fp = join(self.tmpdir, 'tmp.tsv')

        # filter by count
        filter_wf(input_fp, output_fp, min_count=5)
        with open(output_fp, 'r') as f:
            n = len(f.read().splitlines()) - 1
        self.assertEqual(n, 43)

        # filter by percentage
        filter_wf(input_fp, output_fp, min_percent=1)
        with open(output_fp, 'r') as f:
            n = len(f.read().splitlines()) - 1
        self.assertEqual(n, 30)
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


if __name__ == '__main__':
    main()
