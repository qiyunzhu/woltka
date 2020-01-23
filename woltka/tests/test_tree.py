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
from tempfile import mkdtemp, mkstemp
from io import StringIO

from woltka.tree import (
    read_names, read_nodes, read_lineage)


class TreeTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_read_lineage(self):
        map_ = {'seq1': 'k1;p1;c1;o1',
                'seq2': 'k1;p2;c2;o2',
                'seq3': 'k1;p2;c3;',
                'seq4': 'k2;p3;c1;o3',
                'seq5': 'k2;;c4;o4'}
        obs = read_lineage(StringIO('\n'.join([
            f'{k}\t{v}' for k, v in map_.items()])))
        exp = {'k1':          None,
               'k2':          None,
               'k1;p1':       'k1',
               'k1;p2':       'k1',
               'k2;p3':       'k2',
               'k2;':         'k2',
               'k1;p1;c1':    'k1;p1',
               'k1;p2;c2':    'k1;p2',
               'k1;p2;c3':    'k1;p2',
               'k2;p3;c1':    'k2;p3',
               'k2;;c4':      'k2;',
               'k1;p1;c1;o1': 'k1;p1;c1',
               'k1;p2;c2;o2': 'k1;p2;c2',
               'k1;p2;c3;':   'k1;p2;c3',
               'k2;p3;c1;o3': 'k2;p3;c1',
               'k2;;c4;o4':   'k2;;c4'}
        self.assertDictEqual(obs, exp)


if __name__ == '__main__':
    main()
