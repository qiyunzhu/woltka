#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Unit tests for the ordinal system, as well as a demonstration of how the
system works. See docstrings of `test_match_read_gene_dummy`.
"""

from unittest import TestCase, main
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp
from io import StringIO
from functools import partial
from collections import defaultdict

import numpy as np

from woltka.file import openzip
from woltka.align import parse_b6o_file_ext, parse_sam_file_ext
from woltka.ordinal import (
    match_read_gene_dummy, match_read_gene, match_read_gene_naive,
    match_read_gene_quart, ordinal_mapper_dummy, ordinal_mapper, flush_chunk,
    ordinal_mapper_np, flush_chunk_np, load_gene_coords, encode_genes,
    calc_gene_lens, load_gene_lens)


class OrdinalTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_match_read_gene_dummy(self):
        """This test demonstrates how genes and reads are matched in an
        ordinal system. It is for demonstration purpose.

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

        # alignment length map
        lens = {'r{}'.format(i): 20 for i in range(1, 10)}

        # flatten lists
        genes = [x for id_, beg, end in genes for x in (
            (beg,  True,  True, id_),
            (end, False,  True, id_))]
        reads = [x for id_, beg, end in reads for x in (
            (beg,  True, False, id_),
            (end, False, False, id_))]

        # merge genes and reads and sort
        queue = sorted(genes + reads)

        # default (threshold = 80%)
        obs = list(match_read_gene_dummy(queue, lens, th=0.8))
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        self.assertListEqual(obs, exp)

        # threashold = 50%
        obs = list(match_read_gene_dummy(queue, lens, th=0.5))
        exp = [('r1', 'g1'),
               ('r2', 'g1'),
               ('r3', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r7', 'g3'),
               ('r8', 'g3'),
               ('r9', 'g3')]
        self.assertListEqual(obs, exp)

    def test_match_read_gene(self):
        """This test is for the real function.
        """
        # gene table
        genes = [(5, 29),   # 0
                 (33, 61),  # 1
                 (65, 94)]  # 2
        # read map
        reads = [(10, 29),  # 0
                 (16, 35),  # 1
                 (20, 39),  # 2
                 (22, 41),  # 3
                 (30, 49),  # 4
                 (30, 49),  # 5 (identical to 4)
                 (60, 79),  # 6
                 (65, 84),  # 7
                 (82, 95)]  # 8 (shorter than others)

        # read length is uniformly 20, threshold is 80%,
        # so effective length is 20 * 0.8 = 16
        rels = np.full(len(reads), 16, dtype=np.uint16)

        # flatten genes and reads
        genes = np.array([x for idx, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + idx,
            (end << 24) + (3 << 22) + idx)])
        reads = np.array([x for idx, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + idx,
            (end << 24) + (1 << 23) + idx)])
        queue = np.concatenate((genes, reads))
        queue.sort()

        # merge and sort genes and reads
        queue = np.concatenate((genes, reads))
        queue.sort()

        # default protocol
        obs = list(match_read_gene(queue, rels))

        # result: read idx, gene idx
        exp = [(0, 0),
               (4, 1),
               (5, 1),
               (7, 2)]
        self.assertListEqual(obs, exp)

    def test_match_read_gene_naive(self):
        """The naive solution should produce identical result compared to
        the default (ordinal) solution.
        """
        genes = [(5, 29),
                 (33, 61),
                 (65, 94)]
        reads = [(10, 29),
                 (16, 35),
                 (20, 39),
                 (22, 41),
                 (30, 49),
                 (30, 49),
                 (60, 79),
                 (65, 84),
                 (82, 95)]
        # shorten effective length
        rels = np.full(len(reads), 14, dtype=np.uint16)
        genes = np.array([x for idx, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + idx,
            (end << 24) + (3 << 22) + idx)])
        genes.sort()
        reads = [x for idx, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + idx,
            (end << 24) + (1 << 23) + idx)]

        # don't sort, but directly feed both queues
        obs = list(match_read_gene_naive(genes, reads, rels))
        exp = [(0, 0),
               (1, 0),
               (4, 1),
               (5, 1),
               (6, 2),
               (7, 2)]
        self.assertListEqual(obs, exp)

    def test_match_read_gene_quart(self):
        """It should produce identical result compared to the naive method.
        """
        genes = [(5, 29),
                 (33, 61),
                 (65, 94),
                 (61, 76),  # added a small gene within a read
                 (68, 72)]  # added a tiny gene
        reads = [(10, 29),
                 (16, 35),
                 (20, 39),
                 (22, 41),
                 (30, 49),
                 (30, 49),
                 (60, 79),
                 (65, 84),
                 (70, 75),  # added a small read starting in right half
                 (82, 95)]
        rels = np.full(len(reads), 14, dtype=np.uint16)
        genes = np.array([x for idx, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + idx,
            (end << 24) + (3 << 22) + idx)])
        genes.sort()
        reads = [x for idx, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + idx,
            (end << 24) + (1 << 23) + idx)]

        obs = list(match_read_gene_quart(genes, reads, rels))
        exp = [(0, 0),
               (1, 0),
               (4, 1),
               (5, 1),
               (6, 3),  # match to added small gene
               (6, 2),
               (7, 2)]
        self.assertListEqual(obs, exp)

        # a special case with a giant read
        genes = [(1, 5),
                 (6, 7),
                 (7, 8)]
        reads = [(4, 9)]
        rels = np.full(len(reads), 5, dtype=np.uint16)
        genes = np.array([x for idx, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + idx,
            (end << 24) + (3 << 22) + idx)])
        genes.sort()
        reads = np.array([x for idx, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + idx,
            (end << 24) + (1 << 23) + idx)])
        obs = list(match_read_gene_quart(genes, reads, rels))
        exp = []
        self.assertListEqual(obs, exp)

    def test_ordinal_mapper_dummy(self):
        # b6o (BLAST, DIAMOND, BURST, etc.)
        b6o = (
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270')
        obs = ordinal_mapper_dummy(b6o, parse_b6o_file_ext)
        self.assertListEqual(obs[0], ['S1/1', 'S1/2'])
        self.assertDictEqual(obs[1], {'NC_123456': {0: 100, 1: 98}})
        self.assertDictEqual(obs[2], {'NC_123456': [
            (225, True, False, 0), (324, False, False, 0),
            (608, True, False, 1), (708, False, False, 1)]})

        # sam (BWA, Bowtie2, Minimap2 etc.)
        sam = iter((
            # SAM header to be ignored
            '@HD	VN:1.0	SO:unsorted',
            # normal, fully-aligned, forward strand
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            # shortened, reverse strand
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            # not perfectly aligned, unpaired
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            # unaligned
            'S2	16	*	0	0	*	*	0	0	*	*'))
        obs = ordinal_mapper_dummy(sam, parse_sam_file_ext)
        self.assertListEqual(obs[0], ['S1/1', 'S1/2', 'S2'])
        self.assertDictEqual(obs[1], {
            'NC_123456': {0: 100, 1: 80},
            'NC_789012': {2: 90}})
        self.assertDictEqual(obs[2], {
            'NC_123456': [(26,  True, False, 0), (125, False, False, 0),
                          (151, True, False, 1), (230, False, False, 1)],
            'NC_789012': [(186, True, False, 2), (280, False, False, 2)]})

    def test_ordinal_mapper(self):
        # uses the same example as above, with some noises
        coords, idmap, _ = load_gene_coords((
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
        obs = list(ordinal_mapper(aln, coords, idmap))[0]
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        # self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        # self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # # specify format
        # aln.seek(0)
        # obs = list(ordinal_mapper(aln, coords, idmap, fmt='b6o'))[0]
        # self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        # self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify chunk size
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, idmap, n=5))
        self.assertListEqual(list(obs[0][0]), [x[0] for x in exp[:2]])
        self.assertListEqual(list(obs[0][1]), [{x[1]} for x in exp[:2]])
        self.assertListEqual(list(obs[1][0]), [x[0] for x in exp[2:]])
        self.assertListEqual(list(obs[1][1]), [{x[1]} for x in exp[2:]])

        # # add prefix
        # aln.seek(0)
        # obs = list(ordinal_mapper(aln, coords, idmap, prefix=True))[0]
        # self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        # self.assertListEqual(list(obs[1]), [{f'n1_{x[1]}'} for x in exp])

        # # specify threshold
        # aln.seek(0)
        # obs = list(ordinal_mapper(aln, coords, idmap, th=0.5))[0]
        # exp = [('r1', 'g1'),
        #        ('r2', 'g1'),
        #        ('r3', 'g1'),
        #        ('r5', 'g2'),
        #        ('r6', 'g2'),
        #        ('r7', 'g3'),
        #        ('r8', 'g3'),
        #        ('r9', 'g3')]
        # self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        # self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

    def test_flush_chunk(self):
        coords, idmap, _ = load_gene_coords((
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
        idx, rids, rlens = 0, [], []
        locmap = defaultdict(list)
        for idx, row in enumerate(parse_b6o_file_ext(aln)):
            query, subject, _, length, beg, end = row[:6]
            rids.append(query)
            rlens.append(length)
            locmap[subject].extend((
                (beg << 24) + idx, (end << 24) + (1 << 23) + idx))
        rlens = np.array(rlens, dtype=np.uint16)
        obs = flush_chunk(
            idx + 1, locmap, rids, rlens, coords, idmap, 0.8, False)
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

    def test_ordinal_mapper_np(self):
        # should produce the same result as ordinal_mapper
        coords, idmap, _ = load_gene_coords((
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
        obs = list(ordinal_mapper_np(aln, coords, idmap))[0]
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify format
        aln.seek(0)
        obs = list(ordinal_mapper_np(aln, coords, idmap, fmt='b6o'))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify chunk size
        aln.seek(0)
        obs = list(ordinal_mapper_np(aln, coords, idmap, n=5))
        self.assertListEqual(list(obs[0][0]), [x[0] for x in exp[:2]])
        self.assertListEqual(list(obs[0][1]), [{x[1]} for x in exp[:2]])
        self.assertListEqual(list(obs[1][0]), [x[0] for x in exp[2:]])
        self.assertListEqual(list(obs[1][1]), [{x[1]} for x in exp[2:]])

        # add prefix
        aln.seek(0)
        obs = list(ordinal_mapper_np(aln, coords, idmap, prefix=True))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{f'n1_{x[1]}'} for x in exp])

        # specify threshold
        aln.seek(0)
        obs = list(ordinal_mapper_np(aln, coords, idmap, th=0.5))[0]
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

    def test_flush_chunk_np(self):
        coords, idmap, _ = load_gene_coords((
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
        idx = 0
        qrys, subs, lens, begs, ends = [], [], [], [], []
        for idx, row in enumerate(parse_b6o_file_ext(aln)):
            query, subject, _, length, beg, end = row[:6]
            qrys.append(query)
            subs.append(subject)
            lens.append(length)
            begs.append(beg)
            ends.append(end)
        lens = np.array(lens, dtype=np.uint16)
        begs = np.array(begs, dtype=np.int64)
        ends = np.array(ends, dtype=np.int64)
        obs = flush_chunk_np(idx + 1, qrys, subs, lens, begs, ends, coords,
                             idmap, 0.8, False)
        exp = [('r1', 'g1'),
               ('r5', 'g2'),
               ('r6', 'g2'),
               ('r8', 'g3')]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

    def test_load_gene_coords(self):
        # simple case
        tbl = ('>n1',
               'g1	5	29',
               'g2	33	61',
               'g3	65	94')
        obs, idmap, isdup = load_gene_coords(tbl)
        self.assertFalse(isdup)
        self.assertDictEqual(idmap, {'n1': ['g1', 'g2', 'g3']})
        exp = {'n1': np.array([
            (5 << 24) + (1 << 22) + 0,
            (29 << 24) + (3 << 22) + 0,
            (33 << 24) + (1 << 22) + 1,
            (61 << 24) + (3 << 22) + 1,
            (65 << 24) + (1 << 22) + 2,
            (94 << 24) + (3 << 22) + 2])}
        self.assertListEqual(list(obs.keys()), list(exp.keys()))
        np.testing.assert_array_equal(obs['n1'], exp['n1'])

        # NCBI accession
        tbl = ('## GCF_000123456\n',
               '# NC_123456\n',
               '1	5	384\n',
               '2	410	933\n',
               '# NC_789012\n',
               '1	912	638\n',
               '2	529	75\n')
        obs, idmap, isdup = load_gene_coords(tbl, sort=True)
        self.assertTrue(isdup)
        self.assertDictEqual(idmap, {
            'NC_123456': ['1', '2'], 'NC_789012': ['1', '2']})
        exp = {'NC_123456': np.array(
            [(5 << 24) + (1 << 22) + 0,
             (384 << 24) + (3 << 22) + 0,
             (410 << 24) + (1 << 22) + 1,
             (933 << 24) + (3 << 22) + 1]),
               'NC_789012': np.array(
            [(75 << 24) + (1 << 22) + 1,
             (529 << 24) + (3 << 22) + 1,
             (638 << 24) + (1 << 22) + 0,
             (912 << 24) + (3 << 22) + 0])}
        self.assertEqual(obs.keys(), exp.keys())
        for key in obs.keys():
            np.testing.assert_array_equal(obs[key], exp[key])

        # don't sort
        obs = load_gene_coords(tbl, sort=False)[0]['NC_789012']
        exp = np.array([
            (638 << 24) + (1 << 22) + 0, (912 << 24) + (3 << 22) + 0,
            (75 << 24) + (1 << 22) + 1,  (529 << 24) + (3 << 22) + 1])
        np.testing.assert_array_equal(obs, exp)

        # incorrect formats
        # empty file
        msg = 'No coordinate was read from file.'
        with self.assertRaises(ValueError) as ctx:
            load_gene_coords(())
        # only one column
        msg = 'Cannot extract coordinates from line:'
        with self.assertRaises(ValueError) as ctx:
            load_gene_coords(('hello',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello".')
        # only two columns
        with self.assertRaises(ValueError) as ctx:
            load_gene_coords(('hello\t100',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100".')

        # real coords file
        fp = join(self.datdir, 'function', 'coords.txt.xz')
        with openzip(fp) as f:
            obs, idmap, isdup = load_gene_coords(f, sort=True)
        self.assertTrue(isdup)
        self.assertEqual(len(idmap), 107)
        self.assertEqual(len(obs), 107)
        obs_ = obs['G000006745']
        self.assertEqual(len(obs_), 7188)
        self.assertEqual(obs_[0], (372 << 24) + (1 << 22) + 0)
        self.assertEqual(obs_[1], (806 << 24) + (3 << 22) + 0)
        self.assertEqual(obs_[2], (816 << 24) + (1 << 22) + 1)
        self.assertEqual(obs_[3], (2177 << 24) + (3 << 22) + 1)

    def test_encode_genes(self):
        lst = ['5', '384', '410', '933', '912', '638', '529', '75']
        obs = encode_genes(lst)
        exp = np.array([
            (5 << 24) + (1 << 22) + 0, (384 << 24) + (3 << 22) + 0,
            (410 << 24) + (1 << 22) + 1, (933 << 24) + (3 << 22) + 1,
            (638 << 24) + (1 << 22) + 2, (912 << 24) + (3 << 22) + 2,
            (75 << 24) + (1 << 22) + 3,  (529 << 24) + (3 << 22) + 3])
        np.testing.assert_array_equal(obs, exp)

        # coordinate not a number
        with self.assertRaises(ValueError) as ctx:
            encode_genes(['hello', 'there'])
        self.assertEqual(str(ctx.exception), 'Invalid coordinate(s) found.')

    def test_calc_gene_lens(self):
        coords = {'NC_123456': [(5 << 24) + (1 << 22) + 0,
                                (384 << 24) + (3 << 22) + 0,
                                (410 << 24) + (1 << 22) + 1,
                                (933 << 24) + (3 << 22) + 1],
                  'NC_789012': [(75 << 24) + (1 << 22) + 1,
                                (529 << 24) + (3 << 22) + 1,
                                (638 << 24) + (1 << 22) + 0,
                                (912 << 24) + (3 << 22) + 0]}
        idmap = {'NC_123456': ['1', '2'],
                 'NC_789012': ['1', '2']}
        mapper = partial(ordinal_mapper, coords=coords, idmap=idmap,
                         prefix=True)
        obs = calc_gene_lens(mapper)
        exp = {'NC_123456_1': 380,
               'NC_123456_2': 524,
               'NC_789012_2': 455,
               'NC_789012_1': 275}
        self.assertDictEqual(obs, exp)

        idmap = {'NC_123456': ['NP_135792.1', 'NP_246801.2'],
                 'NC_789012': ['NP_258147.1', 'NP_369258.2']}
        mapper = partial(ordinal_mapper, coords=coords, idmap=idmap,
                         prefix=False)
        obs = calc_gene_lens(mapper)
        exp = {'NP_135792.1': 380,
               'NP_246801.2': 524,
               'NP_369258.2': 455,
               'NP_258147.1': 275}
        self.assertDictEqual(obs, exp)

    def test_load_gene_lens(self):
        # simple case
        tbl = ('>n1',
               'g1	5	29',
               'g2	33	61',
               'g3	65	94')
        obs = load_gene_lens(tbl)
        exp = {'g1': 25, 'g2': 29, 'g3': 30}
        self.assertDictEqual(obs, exp)

        # NCBI accession
        tbl = ('## GCF_000123456\n',
               '# NC_123456\n',
               '1	5	384\n',
               '2	410	933\n',
               '# NC_789012\n',
               '1	912	638\n',
               '2	529	75\n')
        obs = load_gene_lens(tbl)
        exp = {'NC_123456_1': 380, 'NC_123456_2': 524,
               'NC_789012_1': 275, 'NC_789012_2': 455}
        self.assertDictEqual(obs, exp)

        # incorrect formats
        # only one column
        msg = 'Cannot calculate gene length from line:'
        with self.assertRaises(ValueError) as ctx:
            load_gene_lens(('hello',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello".')
        # only two columns
        with self.assertRaises(ValueError) as ctx:
            load_gene_lens(('hello\t100',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100".')
        # three columns but 3rd is string
        with self.assertRaises(ValueError) as ctx:
            load_gene_lens(('hello\t100\tthere',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100\tthere".')

        # real coords file
        fp = join(self.datdir, 'function', 'coords.txt.xz')
        with openzip(fp) as f:
            obs = load_gene_lens(f)
        self.assertEqual(obs['G000006845_253'], 996)
        self.assertEqual(obs['G000006745_1171'], 549)
        self.assertEqual(obs['G000006925_4349'], 4224)


if __name__ == '__main__':
    main()
