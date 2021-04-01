#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Unit tests for the ordinal system, as well as a demonstration of how the
system works. See docstrings of `test_match_read_gene`.
"""

from unittest import TestCase, main
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp
from io import StringIO

from woltka.file import openzip
from woltka.align import parse_b6o_line, parse_sam_line
from woltka.ordinal import (
    match_read_gene, match_read_gene_pfx, ordinal_mapper, ordinal_parser,
    read_gene_coords, whether_prefix)


class OrdinalTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_match_read_gene(self):
        """This test demonstrates how genes and reads are matched in an
        ordinal system.

        Illustration of gene (>) and read (-) locations on a genome (=):

           reads:             ---------r1---------
                                    ---------r2---------
                                        ---------r3---------
                                          ---------r4---------
                                                  ---------r5---------
                                                  ---------r6---------
          genome:  1 ====>>>>>>>>>>>g1>>>>>>>>>>>>===>>>>>>>>>>>>>>g2>> 50

           reads:             ---------r7---------
                                   ---------r8---------
                                                    ------r9------
          genome: 51 >>>>>>>>>>>===>>>>>>>>>>>>>>g3>>>>>>>>>>>>>======= 100
          --------------------------------------------------

        This function finds matching read/gene pairs given their overlap
        is no less than a fraction of alignment length (default: 80%). In this
        example, the matching pairs are:

          r1 - g1  (whole read within gene)
          r5 - g2  (17 / 20 nt within gene, exceeding threshold)
          r6 - g2  (same as above)
          r8 - g3  (whole read within gene)

        """
        # gene table
        genes = [('g1',  5, 29),
                 ('g2', 33, 61),
                 ('g3', 65, 94)]
        # read map
        reads = [('r1', 10, 29),
                 ('r2', 16, 35),
                 ('r3', 20, 39),
                 ('r4', 22, 41),
                 ('r5', 30, 49),
                 ('r6', 30, 49),  # identical
                 ('r7', 60, 79),
                 ('r8', 65, 84),
                 ('r9', 82, 95)]  # shorter

        # length map
        lens = {'r{}'.format(i): 20 for i in range(1, 10)}

        # flatten lists
        genes = [x for id_, start, end in genes for x in
                 ((start, True, True, id_),
                  (end,  False, True, id_))]
        reads = [x for id_, start, end in reads for x in
                 ((start, True, False, id_),
                  (end,  False, False, id_))]

        queue = sorted(genes + reads)

        # default (threshold = 80%)
        obs = list(match_read_gene(queue, lens, th=0.8))
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        self.assertListEqual(obs, exp)

        # threashold = 50%
        obs = list(match_read_gene(queue, lens, th=0.5))
        exp = [('r1', 'g1'),
               ('r2', 'g1'),
               ('r3', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r7', 'g3'),
               ('r8', 'g3'),
               ('r9', 'g3')]
        self.assertListEqual(obs, exp)

    def test_match_read_gene_pfx(self):
        # same as above but adds a prefix to genes
        queue = [(1,  True,  True,  'g1'),
                 (11, True,  False, 'r1'),
                 (20, False, False, 'r1'),
                 (26, True,  False, 'r2'),
                 (35, False, False, 'r2'),
                 (50, False, True,  'g1')]
        lens = {'r1': 10, 'r2': 10}
        obs = list(match_read_gene_pfx(queue, lens, th=0.8, pfx='test'))
        exp = [('r1', 'test_g1'),
               ('r2', 'test_g1')]
        self.assertListEqual(obs, exp)

    def test_ordinal_mapper(self):
        # uses the same example as above, with some noises
        coords = read_gene_coords((
            '>n1',
            'g1	5	29',
            'g2	33	61',
            'g3	65	94',
            'gx	108	135'))
        aln = StringIO('\n'.join((
            'r1	n1	95	20	0	0	1	20	10	29	1	1',
            'r2	n1	95	20	0	0	1	20	16	35	1	1',
            'r3	n1	95	20	0	0	1	20	20	39	1	1',
            'r4	n1	95	20	0	0	20	1	22	41	1	1',
            'r5	n1	95	20	0	0	20	1	30	49	1	1',
            'rx	nx	95	20	0	0	1	20	1	20	1	1',
            'r6	n1	95	20	0	0	1	20	49	30	1	1',
            'r7	n1	95	20	0	0	25	6	79	60	1	1',
            'r8	n1	95	20	0	0	1	20	84	65	1	1',
            'r9	n1	95	20	0	0	1	20	95	82	1	1',
            'rx	nx	95	0	0	0	0	0	0	0	1	1',
            '# end of file')))
        obs = list(ordinal_mapper(aln, coords))[0]
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify format
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, fmt='b6o'))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify chunk size
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, n=5))
        self.assertListEqual(list(obs[0][0]), [x[0] for x in exp[:2]])
        self.assertListEqual(list(obs[0][1]), [{x[1]} for x in exp[:2]])
        self.assertListEqual(list(obs[1][0]), [x[0] for x in exp[2:]])
        self.assertListEqual(list(obs[1][1]), [{x[1]} for x in exp[2:]])

        # add prefix
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, prefix=True))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{f'n1_{x[1]}'} for x in exp])

        # specify threshold
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, th=0.5))[0]
        exp = [('r1', 'g1'),
               ('r2', 'g1'),
               ('r3', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r7', 'g3'),
               ('r8', 'g3'),
               ('r9', 'g3')]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

    def test_ordinal_parser(self):

        # b6o (BLAST, DIAMOND, BURST, etc.)
        b6o = (
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270')
        parser = parse_b6o_line
        obs = ordinal_parser(b6o, parser)
        self.assertListEqual(obs[0], ['S1/1', 'S1/2'])
        self.assertDictEqual(obs[1], {'NC_123456': {0: 100, 1: 98}})
        self.assertDictEqual(obs[2], {'NC_123456': [
            (225, True, False, 0), (324, False, False, 0),
            (608, True, False, 1), (708, False, False, 1)]})

        # sam (BWA, Bowtie2, Minimap2 etc.)
        sam = (
            # SAM header to be ignored
            '@HD	VN:1.0	SO:unsorted',
            # normal, fully-aligned, forward strand
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            # shortened, reverse strand
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            # not perfectly aligned, unpaired
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            # unaligned
            'S2	16	*	0	0	*	*	0	0	*	*')
        parser = parse_sam_line
        obs = ordinal_parser(sam, parser)
        self.assertListEqual(obs[0], ['S1/1', 'S1/2', 'S2'])
        self.assertDictEqual(obs[1], {
            'NC_123456': {0: 100, 1: 80},
            'NC_789012': {2: 90}})
        self.assertDictEqual(obs[2], {
            'NC_123456': [(26,  True, False, 0), (125, False, False, 0),
                          (151, True, False, 1), (230, False, False, 1)],
            'NC_789012': [(186, True, False, 2), (280, False, False, 2)]})

    def test_read_gene_coords(self):
        # simple case
        tbl = ('## GCF_000123456\n',
               '# NC_123456\n',
               '1	5	384\n',
               '2	410	933\n',
               '# NC_789012\n',
               '1	912	638\n',
               '2	529	75\n')
        obs = read_gene_coords(tbl, sort=True)
        exp = {'NC_123456': [
            (5,   True, True, '1'), (384, False, True, '1'),
            (410, True, True, '2'), (933, False, True, '2')],
               'NC_789012': [
            (75,  True, True, '2'), (529, False, True, '2'),
            (638, True, True, '1'), (912, False, True, '1')]}
        self.assertDictEqual(obs, exp)

        # don't sort
        obs = read_gene_coords(tbl, sort=False)['NC_789012']
        exp = [(638, True, True, '1'), (912, False, True, '1'),
               (75,  True, True, '2'), (529, False, True, '2')]
        self.assertListEqual(obs, exp)

        # incorrect formats
        # only one column
        msg = 'Cannot extract coordinates from line:'
        with self.assertRaises(ValueError) as ctx:
            read_gene_coords(('hello',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello".')
        # only two columns
        with self.assertRaises(ValueError) as ctx:
            read_gene_coords(('hello\t100',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100".')
        # three columns but 3rd is string
        with self.assertRaises(ValueError) as ctx:
            read_gene_coords(('hello\t100\tthere',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100\tthere".')

        # real coords file
        fp = join(self.datdir, 'function', 'coords.txt.xz')
        with openzip(fp) as f:
            obs = read_gene_coords(f, sort=True)
        self.assertEqual(len(obs), 107)
        obs_ = obs['G000006745']
        self.assertEqual(len(obs_), 7188)
        self.assertTupleEqual(obs_[0], (372,  True,  True, '1'))
        self.assertTupleEqual(obs_[1], (806,  False, True, '1'))
        self.assertTupleEqual(obs_[2], (816,  True,  True, '2'))
        self.assertTupleEqual(obs_[3], (2177, False, True, '2'))

    def test_whether_prefix(self):
        # gene Ids are indices
        coords = {'NC_123456': [
            (5,   True, True, '1'), (384, False, True, '1'),
            (410, True, True, '2'), (933, False, True, '2')],
                  'NC_789012': [
            (75,  True, True, '2'), (529, False, True, '2'),
            (638, True, True, '1'), (912, False, True, '1')]}
        self.assertTrue(whether_prefix(coords))

        # gene Ids are unique accessions
        coords = {'NC_123456': [
            (5,   True, True,  'NP_135792.1'),
            (384, False, True, 'NP_135792.1'),
            (410, True, True,  'NP_246801.2'),
            (933, False, True, 'NP_246801.2')],
                  'NC_789012': [
            (75,  True, True,  'NP_258147.1'),
            (529, False, True, 'NP_258147.1'),
            (638, True, True,  'NP_369258.2'),
            (912, False, True, 'NP_369258.2')]}
        self.assertFalse(whether_prefix(coords))


if __name__ == '__main__':
    main()
