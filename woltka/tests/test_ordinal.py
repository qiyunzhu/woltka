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

from woltka.file import openzip
from woltka.align import parse_b6o_file_ex, parse_sam_file_ex
from woltka.ordinal import (
    match_read_gene_dummy, match_read_gene, match_read_gene_naive,
    match_read_gene_quart, ordinal_mapper, flush_chunk,
    load_gene_coords, calc_gene_lens, load_gene_lens)


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
        genes = [('g1',  4, 29),
                 ('g2', 32, 61),
                 ('g3', 64, 94)]
        # read map
        reads = [('r1',  9, 29),
                 ('r2', 15, 35),
                 ('r3', 19, 39),
                 ('r4', 21, 41),
                 ('r5', 29, 49),
                 ('r6', 29, 49),  # identical
                 ('r7', 59, 79),
                 ('r8', 64, 84),
                 ('r9', 81, 95)]  # shorter

        # alignment length map
        lens = {'r{}'.format(i): 20 for i in range(1, 10)}

        # flatten lists
        genes = [x for id_, beg, end in genes for x in (
            (beg, False,  True, id_),
            (end,  True,  True, id_))]
        reads = [x for id_, beg, end in reads for x in (
            (beg, False, False, id_),
            (end,  True, False, id_))]

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

        # read length is uniformly 20, threshold is 80%, so effective
        # alignment length is 20 * 0.8 = 16
        rels = [16] * len(reads)
        genes = [x for i, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + i,
            (end << 24) + (3 << 22) + i)]
        reads = [x for i, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + i,
            (end << 24) + (1 << 23) + i)]
        queue = sorted(genes + reads)

        # default
        obs = list(match_read_gene(queue, rels))
        exp = [(0, 0),
               (4, 1),
               (5, 1),
               (7, 2)]
        self.assertListEqual(obs, exp)

    def test_match_read_gene_naive(self):
        """The naive solution should produce identical result compared to
        the default (ordinal) solution.
        """
        genes = [(4, 29),
                 (32, 61),
                 (64, 94)]
        reads = [(9, 29),
                 (15, 35),
                 (19, 39),
                 (21, 41),
                 (29, 49),
                 (29, 49),
                 (59, 79),
                 (64, 84),
                 (81, 95)]
        # shorten effective length
        rels = [14] * len(reads)
        genes = [x for i, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + i,
            (end << 24) + (3 << 22) + i)]
        reads = [x for i, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + i,
            (end << 24) + (1 << 23) + i)]

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
        """The naive solution should produce identical result compared to
        the default (ordinal) solution.
        """
        genes = [(4, 29),
                 (32, 61),
                 (64, 94),
                 (60, 76),  # added a small gene within a read
                 (66, 72)]  # added a tiny gene
        reads = [(9, 29),
                 (15, 35),
                 (19, 39),
                 (21, 41),
                 (29, 49),
                 (29, 49),
                 (59, 79),
                 (64, 84),
                 (69, 75),  # added a small read starting in right half
                 (81, 95)]
        rels = [14] * len(reads)
        genes = [x for i, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + i,
            (end << 24) + (3 << 22) + i)]
        genes.sort()
        reads = [x for i, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + i,
            (end << 24) + (1 << 23) + i)]

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
        genes = [(0, 5),
                 (5, 7),
                 (6, 8)]
        reads = [(3, 9)]
        rels = [5]
        genes = [x for i, (beg, end) in enumerate(genes) for x in (
            (beg << 24) + (1 << 22) + i,
            (end << 24) + (3 << 22) + i)]
        genes.sort()
        reads = [x for i, (beg, end) in enumerate(reads) for x in (
            (beg << 24) + i,
            (end << 24) + (1 << 23) + i)]
        obs = list(match_read_gene_quart(genes, reads, rels))
        exp = []
        self.assertListEqual(obs, exp)

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
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify format
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, idmap, fmt='b6o'))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{x[1]} for x in exp])

        # specify chunk size
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, idmap, n=5))
        self.assertListEqual(list(obs[0][0]), [x[0] for x in exp[:2]])
        self.assertListEqual(list(obs[0][1]), [{x[1]} for x in exp[:2]])
        self.assertListEqual(list(obs[1][0]), [x[0] for x in exp[2:]])
        self.assertListEqual(list(obs[1][1]), [{x[1]} for x in exp[2:]])

        # add prefix
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, idmap, prefix=True))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [{f'n1_{x[1]}'} for x in exp])

        # specify threshold
        aln.seek(0)
        obs = list(ordinal_mapper(aln, coords, idmap, th=0.5))[0]
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
        idx, rids, rlens, locmap = 0, [], [], {}
        for query, records in parse_b6o_file_ex(aln):
            for subject, _, length, beg, end in records:
                rids.append(query)
                rlens.append(length)
                locmap.setdefault(subject, []).extend((
                    (beg << 24) + idx, (end << 24) + (1 << 23) + idx))
                idx += 1
        obs = flush_chunk(locmap, rids, rlens, coords, idmap, 0.8, False)
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
        exp = {'n1': [
            (4 << 24) + (1 << 22) + 0,
            (29 << 24) + (3 << 22) + 0,
            (32 << 24) + (1 << 22) + 1,
            (61 << 24) + (3 << 22) + 1,
            (64 << 24) + (1 << 22) + 2,
            (94 << 24) + (3 << 22) + 2]}
        self.assertDictEqual(obs, exp)

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
        exp = {'NC_123456': [
            (4 << 24) + (1 << 22) + 0,
            (384 << 24) + (3 << 22) + 0,
            (409 << 24) + (1 << 22) + 1,
            (933 << 24) + (3 << 22) + 1],
               'NC_789012': [
            (74 << 24) + (1 << 22) + 1,
            (529 << 24) + (3 << 22) + 1,
            (637 << 24) + (1 << 22) + 0,
            (912 << 24) + (3 << 22) + 0]}
        self.assertDictEqual(obs, exp)

        # don't sort
        obs = load_gene_coords(tbl, sort=False)[0]['NC_789012']
        exp = [(637 << 24) + (1 << 22) + 0, (912 << 24) + (3 << 22) + 0,
               (74 << 24) + (1 << 22) + 1,  (529 << 24) + (3 << 22) + 1]
        self.assertListEqual(obs, exp)

        # incorrect formats
        # only one column
        msg = 'Cannot extract coordinates from line:'
        with self.assertRaises(ValueError) as ctx:
            load_gene_coords(('hello',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello".')
        # only two columns
        with self.assertRaises(ValueError) as ctx:
            load_gene_coords(('hello\t100',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100".')
        # three columns but 3rd is string
        with self.assertRaises(ValueError) as ctx:
            load_gene_coords(('hello\t100\tthere',))
        self.assertEqual(str(ctx.exception), f'{msg} "hello\t100\tthere".')

        # real coords file
        fp = join(self.datdir, 'function', 'coords.txt.xz')
        with openzip(fp) as f:
            obs, idmap, isdup = load_gene_coords(f, sort=True)
        self.assertTrue(isdup)
        self.assertEqual(len(idmap), 107)
        self.assertEqual(len(obs), 107)
        obs_ = obs['G000006745']
        self.assertEqual(len(obs_), 7188)
        self.assertEqual(obs_[0], (371 << 24) + (1 << 22) + 0)
        self.assertEqual(obs_[1], (806 << 24) + (3 << 22) + 0)
        self.assertEqual(obs_[2], (815 << 24) + (1 << 22) + 1)
        self.assertEqual(obs_[3], (2177 << 24) + (3 << 22) + 1)

    def test_calc_gene_lens(self):
        coords = {'NC_123456': [(4 << 24) + (1 << 22) + 0,
                                (384 << 24) + (3 << 22) + 0,
                                (409 << 24) + (1 << 22) + 1,
                                (933 << 24) + (3 << 22) + 1],
                  'NC_789012': [(74 << 24) + (1 << 22) + 1,
                                (529 << 24) + (3 << 22) + 1,
                                (637 << 24) + (1 << 22) + 0,
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
