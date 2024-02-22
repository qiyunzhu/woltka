#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp
from io import StringIO

from woltka.file import openzip
from woltka.align import (
    plain_mapper, range_mapper, iter_align, infer_align_format, assign_parser,
    parse_sam_file, parse_sam_file_ex, parse_sam_file_ft, parse_sam_file_ex_ft,
    cigar_to_lens, cigar_to_lens_ord, parse_sam_file_pd,
    parse_map_file, parse_map_file_ft,
    parse_b6o_file, parse_b6o_file_ex, parse_b6o_file_ft, parse_b6o_file_ex_ft,
    parse_paf_file, parse_paf_file_ex, parse_paf_file_ft, parse_paf_file_ex_ft)


class AlignTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_plain_mapper(self):

        def _res2lst(res):
            return tuple(tuple(list(x) for x in y) for y in res)

        # simple case
        aln = StringIO('\n'.join((
            'R1	G1',
            'R2	G1',
            'R2	G2',
            'R3	G1',
            'R3	G3',
            'R4	G4',
            'R5	G5')))
        obs = plain_mapper(aln)
        exp = ((['R1', 'R2', 'R3', 'R4', 'R5'],
                [{'G1'}, {'G1', 'G2'}, {'G1', 'G3'}, {'G4'}, {'G5'}]),)
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 1
        aln.seek(0)
        obs = plain_mapper(aln, n=1)
        exp = ((['R1'], [{'G1'}]),
               (['R2'], [{'G1', 'G2'}]),
               (['R3'], [{'G1', 'G3'}]),
               (['R4'], [{'G4'}]),
               (['R5'], [{'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 2
        aln.seek(0)
        obs = plain_mapper(aln, n=2)
        exp = ((['R1', 'R2'], [{'G1'}, {'G1', 'G2'}]),
               (['R3', 'R4'], [{'G1', 'G3'}, {'G4'}]),
               (['R5'], [{'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 3
        aln.seek(0)
        obs = plain_mapper(aln, n=3)
        exp = ((['R1', 'R2', 'R3'], [{'G1'}, {'G1', 'G2'}, {'G1', 'G3'}]),
               (['R4', 'R5'], [{'G4'}, {'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 4
        aln.seek(0)
        obs = plain_mapper(aln, n=4)
        exp = ((['R1', 'R2', 'R3', 'R4'],
                [{'G1'}, {'G1', 'G2'}, {'G1', 'G3'}, {'G4'}]),
               (['R5'], [{'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 5
        aln.seek(0)
        obs = plain_mapper(aln, n=5)
        exp = ((['R1', 'R2', 'R3', 'R4', 'R5'],
                [{'G1'}, {'G1', 'G2'}, {'G1', 'G3'}, {'G4'}, {'G5'}]),)
        self.assertTupleEqual(_res2lst(obs), exp)

        # format is given
        aln.seek(0)
        obs = plain_mapper(aln, fmt='map', n=5)
        self.assertTupleEqual(_res2lst(obs), exp)

        # exclude a subject
        aln.seek(0)
        obs = plain_mapper(aln, excl={'G1'})
        exp = ((['R4', 'R5'], [{'G4'}, {'G5'}]),)
        self.assertTupleEqual(_res2lst(obs), exp)

        # empty alignment
        with self.assertRaises(ValueError) as ctx:
            list(plain_mapper(StringIO()))
        self.assertEqual(str(ctx.exception), (
            'Alignment file is empty or unreadable.'))

        # bad alignment
        with self.assertRaises(ValueError) as ctx:
            list(plain_mapper(StringIO('Hi there!')))
        self.assertEqual(str(ctx.exception), (
            'Cannot determine alignment file format.'))

    def test_range_mapper(self):

        def _res2lst(res):
            return tuple(tuple(list(x) for x in y) for y in res)

        aln = StringIO('\n'.join((
            'R1	G1	95	20	0	0	1	20	10	29	1	1',
            'R2	G1	95	20	0	0	1	20	16	35	1	1',
            'R3	G2	95	20	0	0	1	20	39	21	1	1',
            'R3	G3	95	20	0	0	1	20	88	70	1	1',
            'R4	G2	95	20	0	0	20	1	41	22	1	1',
            'R5	G3	95	20	0	0	20	1	30	49	1	1',
            'R5	G3	95	20	0	0	20	1	50	69	1	1',
            'Rx	Gx	95	20	0	0	1	20	0	0	1	1',
            '# this is not an alignment')))
        obs = _res2lst(range_mapper(aln))[0]
        exp = [('R1', {'G1': [10, 29]}),
               ('R2', {'G1': [16, 35]}),
               ('R3', {'G2': [21, 39], 'G3': [70, 88]}),
               ('R4', {'G2': [22, 41]}),
               ('R5', {'G3': [30, 49, 50, 69]})]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [x[1] for x in exp])

        # chunk of 3
        aln.seek(0)
        obs = _res2lst(range_mapper(aln, n=3))
        self.assertListEqual(list(obs[0][0]), [x[0] for x in exp[:3]])
        self.assertListEqual(list(obs[0][1]), [x[1] for x in exp[:3]])
        self.assertListEqual(list(obs[1][0]), [x[0] for x in exp[3:]])
        self.assertListEqual(list(obs[1][1]), [x[1] for x in exp[3:]])

        # specify format
        aln.seek(0)
        obs = _res2lst(range_mapper(aln, fmt='b6o'))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [x[1] for x in exp])

        # exclude a subject
        aln.seek(0)
        obs = _res2lst(range_mapper(aln, excl={'G1'}))[0]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp[2:]])
        self.assertListEqual(list(obs[1]), [x[1] for x in exp[2:]])

    def test_iter_align(self):
        line = 'S1/1	NC_123456'
        obs = list(iter_align(StringIO(line)))
        self.assertTupleEqual(obs[0], ('S1/1', {'NC_123456'}))
        obs = list(iter_align(StringIO(line), fmt='map'))
        self.assertTupleEqual(obs[0], ('S1/1', {'NC_123456'}))
        obs = list(iter_align(StringIO(line), excl={'abc'}))
        self.assertTupleEqual(obs[0], ('S1/1', {'NC_123456'}))
        obs = list(iter_align(StringIO(line), extr=True))
        self.assertTupleEqual(obs[0], ('S1/1', {'NC_123456'}))

    def test_infer_align_format(self):
        # simple cases
        # map
        line = 'S1/1	NC_123456'
        obs = infer_align_format(StringIO(line))
        self.assertEqual(obs[0], 'map')
        self.assertListEqual(obs[1], [line])

        # b6o
        line = 'S1/1	NC_123456	100	100	0	0	1	100	25	124	0.1	100'
        self.assertEqual(infer_align_format(StringIO(line))[0], 'b6o')

        # paf
        line = ('S1/1	100	0	75	+	NC_123456	1000	100	175	70	75	20'
                '	tp:A:P	cm:i:10	s1:i:100	s2:i:0	rl:i:0')
        self.assertEqual(infer_align_format(StringIO(line))[0], 'paf')

        # sam
        line = 'S1	77	NC_123456	26	0	100M	*	0	0	*	*'
        self.assertEqual(infer_align_format(StringIO(line))[0], 'sam')

        # sam header
        line = '@HD	VN:1.0	SO:unsorted'
        self.assertEqual(infer_align_format(StringIO(line))[0], 'sam')

        # empty file
        with self.assertRaises(ValueError) as ctx:
            infer_align_format(StringIO())
        self.assertEqual(str(ctx.exception), (
            'Alignment file is empty or unreadable.'))

        # invalid sam
        line = 'S1	*	*	*	*	*	*	0	0	*	*'
        with self.assertRaises(ValueError) as ctx:
            infer_align_format(StringIO(line))
        self.assertEqual(str(ctx.exception), (
            'Cannot determine alignment file format.'))

        # cannot determine
        line = 'Hi there!'
        with self.assertRaises(ValueError) as ctx:
            infer_align_format(StringIO(line))
        self.assertEqual(str(ctx.exception), (
            'Cannot determine alignment file format.'))

        # real files
        # Bowtie2 (sam)
        with openzip(join(self.datdir, 'align', 'bowtie2', 'S01.sam.xz')) as f:
            self.assertEqual(infer_align_format(f)[0], 'sam')

        # BURST (b6o)
        with openzip(join(self.datdir, 'align', 'burst', 'S01.b6.bz2')) as f:
            self.assertEqual(infer_align_format(f)[0], 'b6o')

    def test_assign_parser(self):
        self.assertEqual(assign_parser('map'), parse_map_file)
        self.assertEqual(assign_parser('sam'), parse_sam_file)
        self.assertEqual(assign_parser('paf'), parse_paf_file)
        self.assertEqual(assign_parser('b6o'), parse_b6o_file)
        self.assertEqual(assign_parser('map', extr=True), parse_map_file)
        self.assertEqual(assign_parser('sam', extr=True), parse_sam_file_ex)
        self.assertEqual(assign_parser('paf', extr=True), parse_paf_file_ex)
        self.assertEqual(assign_parser('b6o', extr=True), parse_b6o_file_ex)
        self.assertEqual(assign_parser('map', filt=True), parse_map_file_ft)
        self.assertEqual(assign_parser('sam', filt=True), parse_sam_file_ft)
        self.assertEqual(assign_parser('paf', filt=True), parse_paf_file_ft)
        self.assertEqual(assign_parser('b6o', filt=True), parse_b6o_file_ft)
        self.assertEqual(assign_parser('map', extr=True, filt=True),
                         parse_map_file_ft)
        self.assertEqual(assign_parser('sam', extr=True, filt=True),
                         parse_sam_file_ex_ft)
        self.assertEqual(assign_parser('paf', extr=True, filt=True),
                         parse_paf_file_ex_ft)
        self.assertEqual(assign_parser('b6o', extr=True, filt=True),
                         parse_b6o_file_ex_ft)
        with self.assertRaises(ValueError) as ctx:
            assign_parser('xyz')
        self.assertEqual(str(ctx.exception), (
            'Invalid format code: "xyz".'))

    def test_parse_sam_file(self):
        sam = (
            # 1st line: header
            '@HD	VN:1.0	SO:unsorted',
            # 2nd line: normal, fully-aligned, forward strand
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            # 3rd line: shortened, reverse strand
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            # 4th line: not perfectly aligned, unpaired
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            # 5th line: not aligned
            'S2	16	*	0	0	*	*	0	0	*	*',
            # 6-9th lines: interleaved forward and reverse reads
            'S3	83	NC_123456	452	0	100M	*	0	0	*	*',
            'S3	163	NC_123456	378	0	80M5D15M	*	0	0	*	*',
            'S3	355	NC_345678	133	0	100M	*	0	0	*	*',
            'S3	403	NC_345678	261	0	10M5I85M	*	0	0	*	*'
        )
        obs = list(parse_sam_file(iter(sam)))
        exp = [
            ('S1/1', {'NC_123456'}),
            ('S1/2', {'NC_123456'}),
            ('S2',   {'NC_789012'}),
            ('S3/1', {'NC_123456', 'NC_345678'}),
            ('S3/2', {'NC_123456', 'NC_345678'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)

        # header only
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            '@SQ	SN:G000005825	LN:4249288',
            '@SQ	SN:G000006175	LN:1936387'
        )
        self.assertFalse(list(parse_sam_file(iter(sam))))

        # empty file
        self.assertFalse(list(parse_sam_file(iter(()))))

    def test_parse_sam_file_ex(self):
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            'S2	16	*	0	0	*	*	0	0	*	*',
            'S3	83	NC_123456	452	0	100M	*	0	0	*	*',
            'S3	163	NC_123456	378	0	80M5D15M	*	0	0	*	*',
            'S3	355	NC_345678	133	0	100M	*	0	0	*	*',
            'S3	403	NC_345678	261	0	10M5I85M	*	0	0	*	*'
        )
        obs = list(parse_sam_file_ex(iter(sam)))
        exp = [
            ('S1/1', [('NC_123456', None, 100, 26,  125)]),
            ('S1/2', [('NC_123456', None,  80, 151, 230)]),
            ('S2',   [('NC_789012', None,  90, 186, 280)]),
            ('S3/1', [('NC_123456', None, 100, 452, 551),
                      ('NC_345678', None, 100, 133, 232)]),
            ('S3/2', [('NC_123456', None,  95, 378, 477),
                      ('NC_345678', None,  95, 261, 355)])
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        self.assertFalse(list(parse_sam_file_ex(iter(()))))

    def test_parse_sam_file_ft(self):
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	G1	80	0	50M	*	0	0	*	*',
            'S1	141	G2	80	0	50M	*	0	0	*	*',
            'S2	0	G2	80	0	50M	*	0	0	*	*',
            'S2	16	G3	80	0	50M	*	0	0	*	*',
            'S2	147	G4	80	0	50M	*	0	0	*	*',
            'S2	99	G3	80	0	50M	*	0	0	*	*',
            'S2	0	*	0	0	*	*	0	0	*	*',
            'S3	83	G3	80	0	50M	*	0	0	*	*',
            'S3	163	G1	80	0	50M	*	0	0	*	*',
            'S3	355	G4	80	0	50M	*	0	0	*	*',
            'S3	403	G5	80	0	50M	*	0	0	*	*',
            'S4	99	G6	80	0	50M	*	0	0	*	*',
            'S4	147	G6	80	0	50M	*	0	0	*	*',
            'S4	256	G6	80	0	50M	*	0	0	*	*'
        )
        obs = list(parse_sam_file_ft(iter(sam), {'G1'}))
        exp = [
            ('S2',   {'G2', 'G3'}),
            ('S2/1', {'G3'}),
            ('S2/2', {'G4'}),
            ('S4',   {'G6'}),
            ('S4/1', {'G6'}),
            ('S4/2', {'G6'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)

        # single record
        sam = ('S1	77	G1	80	0	50M	*	0	0	*	*',)
        self.assertFalse(list(parse_sam_file_ft(iter(sam), {'G1'})))

        # empty file
        self.assertFalse(list(parse_sam_file_ft(iter(()), {'G1'})))

    def test_parse_sam_file_ex_ft(self):
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	G1	80	0	50M	*	0	0	*	*',
            'S1	141	G2	80	0	50M	*	0	0	*	*',
            'S2	0	G2	80	0	50M	*	0	0	*	*',
            'S2	16	G3	80	0	50M	*	0	0	*	*',
            'S2	147	G4	80	0	50M	*	0	0	*	*',
            'S2	99	G3	80	0	50M	*	0	0	*	*',
            'S2	0	*	0	0	*	*	0	0	*	*',
            'S3	83	G3	80	0	50M	*	0	0	*	*',
            'S3	163	G1	80	0	50M	*	0	0	*	*',
            'S3	355	G4	80	0	50M	*	0	0	*	*',
            'S3	403	G5	80	0	50M	*	0	0	*	*',
            'S4	99	G6	80	0	50M	*	0	0	*	*',
            'S4	147	G6	80	0	50M	*	0	0	*	*',
            'S4	256	G6	80	0	50M	*	0	0	*	*'
        )
        obs = list(parse_sam_file_ex_ft(iter(sam), {'G1'}))
        exp = [
            ('S2',   [('G2', None, 50, 80, 129),
                      ('G3', None, 50, 80, 129)]),
            ('S2/1', [('G3', None, 50, 80, 129)]),
            ('S2/2', [('G4', None, 50, 80, 129)]),
            ('S4',   [('G6', None, 50, 80, 129)]),
            ('S4/1', [('G6', None, 50, 80, 129)]),
            ('S4/2', [('G6', None, 50, 80, 129)])
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        sam = ('S1	77	G1	80	0	50M	*	0	0	*	*',)
        self.assertFalse(list(parse_sam_file_ex_ft(iter(sam), {'G1'})))
        self.assertFalse(list(parse_sam_file_ex_ft(iter(()), {'G1'})))

    def test_cigar_to_lens(self):
        self.assertTupleEqual(cigar_to_lens('150M'), (150, 150))
        self.assertTupleEqual(cigar_to_lens('3M1I3M1D5M'), (11, 12))

    def test_cigar_to_lens_ord(self):
        self.assertTupleEqual(cigar_to_lens_ord('150M'), (150, 150))
        self.assertTupleEqual(cigar_to_lens_ord('3M1I3M1D5M'), (11, 12))

    def test_parse_sam_file_pd(self):
        self.assertIsNone(parse_sam_file_pd([]))

    def test_parse_map_file(self):
        aln = (
            'R1	G1',          # 1st line (normal)
            'R1	G2',          # 2nd line (normal, same query)
            'R2	G3	# note',  # 3rd line (with additional column)
            '',               # 4th line (empty)
            'R2	G4',          # 5th line (same query)
            '# end'           # 6th line (not tab-delimited)
        )
        obs = list(parse_map_file(iter(aln)))
        exp = [
            ('R1', {'G1', 'G2'}),
            ('R2', {'G3', 'G4'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        self.assertFalse(list(parse_map_file(iter(('hello')))))
        self.assertFalse(list(parse_map_file(iter(()))))

    def test_parse_map_file_ft(self):
        aln = (
            'R1	G1',
            'R1	G2',
            'R2	G3	# note',
            '',
            'R2	G4',
            '# end'
        )

        # exclude nothing
        obs = list(parse_map_file_ft(iter(aln), set()))
        exp = [
            ('R1', {'G1', 'G2'}),
            ('R2', {'G3', 'G4'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)

        # exclude 1st subject of 1st query, and 2nd (i.e., non-1st) subject of
        # 2nd (i.e., non-1st) query
        obs = list(parse_map_file_ft(iter(aln), {'G1', 'G4'}))
        self.assertFalse(obs)

        # exclude 2nd subject of 1st query, and 1st subject of 2nd query
        obs = list(parse_map_file_ft(iter(aln), {'G2', 'G3'}))
        self.assertFalse(obs)

        # no valid mapping
        obs = list(parse_map_file_ft(iter(('hello')), {'X'}))
        self.assertFalse(obs)

    def test_parse_b6o_file(self):
        b6o = (
            '# BLAST result:',
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270',
            'S1/2	NC_789012	80	95	5	2	3	97	337	425	5.6e-15	206',
            '# end'
        )
        obs = list(parse_b6o_file(iter(b6o)))
        exp = [
            ('S1/1', {'NC_123456'}),
            ('S1/2', {'NC_123456', 'NC_789012'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        self.assertFalse(list(parse_b6o_file(iter(('hello')))))
        self.assertFalse(list(parse_b6o_file(iter(()))))

    def test_parse_b6o_file_ex(self):
        b6o = (
            '# BLAST result:',
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270',
            'S1/2	NC_789012	80	95	5	2	3	97	337	425	5.6e-15	206',
            '# end'
        )
        obs = list(parse_b6o_file_ex(iter(b6o)))
        self.assertEqual(len(obs), 2)
        exp = [
            ('S1/1', [('NC_123456', 345, 100, 225, 324)]),
            ('S1/2', [('NC_123456', 270, 98, 608, 708),
                      ('NC_789012', 206, 95, 337, 425)])
        ]
        self.assertTupleEqual(obs[0], exp[0])
        self.assertTupleEqual(obs[1], exp[1])

        obs = list(parse_b6o_file_ex(iter(('Hi there!'))))
        self.assertEqual(len(obs), 0)

    def test_parse_b6o_file_ft(self):
        b6o = (
            'S1	G1	100	100	0	0	1	100	1	100	1.1e-9	900',
            'S2	G1	90	90	1	1	5	95	101	200	2.2e-8	800',
            'S2	G2	80	80	2	2	10	90	201	300	3.3e-7	700',
            'S3	G3	70	70	3	3	15	85	301	400	4.4e-6	600',
            'S3	G1	60	60	4	4	20	80	401	500	5.5e-5	500',
            'S4	G2	50	50	5	5	25	75	501	600	6.6e-4	400',
            'S4	G3	40	40	6	6	30	70	601	700	7.7e-5	300',
            '# end'
        )
        obs = list(parse_b6o_file_ft(iter(b6o), set()))
        exp = [
            ('S1', {'G1'}),
            ('S2', {'G1', 'G2'}),
            ('S3', {'G3', 'G1'}),
            ('S4', {'G2', 'G3'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        obs = list(parse_b6o_file_ft(iter(b6o), {'G1'}))
        self.assertEqual(len(obs), 1)
        self.assertTupleEqual(obs[0], ('S4', {'G2', 'G3'}))
        self.assertFalse(list(parse_b6o_file_ft(iter(('hello')), {'G1'})))
        self.assertFalse(list(parse_b6o_file_ft(iter(()), {'G1'})))

    def test_parse_b6o_file_ex_ft(self):
        b6o = (
            'S1	G1	100	100	0	0	1	100	1	100	1.1e-9	900',
            'S2	G1	90	90	1	1	5	95	101	200	2.2e-8	800',
            'S2	G2	80	80	2	2	10	90	201	300	3.3e-7	700',
            'S3	G3	70	70	3	3	15	85	301	400	4.4e-6	600',
            'S3	G1	60	60	4	4	20	80	401	500	5.5e-5	500',
            'S4	G2	50	50	5	5	25	75	501	600	6.6e-4	400',
            'S4	G3	40	40	6	6	30	70	601	700	7.7e-5	300',
            '# end'
        )
        obs = list(parse_b6o_file_ex_ft(iter(b6o), set()))
        exp = [
            ('S1', [('G1', 900, 100,   1, 100)]),
            ('S2', [('G1', 800,  90, 101, 200),
                    ('G2', 700,  80, 201, 300)]),
            ('S3', [('G3', 600,  70, 301, 400),
                    ('G1', 500,  60, 401, 500)]),
            ('S4', [('G2', 400,  50, 501, 600),
                    ('G3', 300,  40, 601, 700)])
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        obs = list(parse_b6o_file_ex_ft(iter(b6o), {'G1'}))
        self.assertEqual(len(obs), 1)
        self.assertTupleEqual(obs[0], exp[3])
        self.assertFalse(list(parse_b6o_file_ex_ft(iter(('hello')), {'G1'})))
        self.assertFalse(list(parse_b6o_file_ex_ft(iter(()), {'G1'})))

    def test_parse_paf_file(self):
        paf = (
            '# Minimap2 result:',
            'S1/1	150	25	100	+	NC_123456	1000000	1025	1100	70	75	20'
            '	tp:A:P	cm:i:10	s1:i:100	s2:i:0	rl:i:0',
            'S1/2	150	25	125	-	NC_123456	1000000	1200	1300	90	100	6'
            '	tp:A:P	cm:i:15	s1:i:150	s2:i:120	rl:i:0',
            'S2/1	250	50	200	-	NC_789012	500000	100000	100150	100	150	65'
            '	tp:A:P	cm:i:24	s1:i:214	s2:i:193	rl:i:0',
            'S2/1	250	50	150	+	NC_345678	1500000	125000	125100	75	100	50'
            '	tp:A:P	cm:i:7	s1:i:55	s2:i:0	rl:i:0',
            'S2/2	250	0	250	+	NC_789012	500000	99900	100150	245	250	0'
            '	tp:A:P	cm:i:24	s1:i:282	s2:i:267	rl:i:0',
            '- the end -'
        )
        obs = list(parse_paf_file(iter(paf)))
        exp = [
            ('S1/1', {'NC_123456'}),
            ('S1/2', {'NC_123456'}),
            ('S2/1', {'NC_789012', 'NC_345678'}),
            ('S2/2', {'NC_789012'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        self.assertFalse(list(parse_paf_file(iter(('hello')))))
        self.assertFalse(list(parse_paf_file(iter(()))))

    def test_parse_paf_file_ex(self):
        paf = (
            '# Minimap2 result:',
            'S1/1	150	25	100	+	NC_123456	1000000	1025	1100	70	75	20'
            '	tp:A:P	cm:i:10	s1:i:100	s2:i:0	rl:i:0',
            'S1/2	150	25	125	-	NC_123456	1000000	1200	1300	90	100	6'
            '	tp:A:P	cm:i:15	s1:i:150	s2:i:120	rl:i:0',
            'S2/1	250	50	200	-	NC_789012	500000	100000	100150	100	150	65'
            '	tp:A:P	cm:i:24	s1:i:214	s2:i:193	rl:i:0',
            'S2/1	250	50	150	+	NC_345678	1500000	125000	125100	75	100	50'
            '	tp:A:P	cm:i:7	s1:i:55	s2:i:0	rl:i:0',
            'S2/2	250	0	250	+	NC_789012	500000	99900	100150	245	250	0'
            '	tp:A:P	cm:i:24	s1:i:282	s2:i:267	rl:i:0',
            '- the end -'
        )
        obs = list(parse_paf_file_ex(iter(paf)))
        exp = [
            ('S1/1', [('NC_123456', 20, 75, 1026, 1100)]),
            ('S1/2', [('NC_123456', 6, 100, 1201, 1300)]),
            ('S2/1', [('NC_789012', 65, 150, 100001, 100150),
                      ('NC_345678', 50, 100, 125001, 125100)]),
            ('S2/2', [('NC_789012', 0, 250, 99901, 100150)])
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        self.assertFalse(list(parse_paf_file_ex(iter(('hello')))))
        self.assertFalse(list(parse_paf_file_ex(iter(()))))

    def test_parse_paf_file_ft(self):
        paf = (
            'S1	100	0	100	+	G1	1000	0	100	100	100	200',
            'S2	100	5	95	-	G1	1000	100	200	90	100	175',
            'S2	100	10	90	+	G2	2000	200	300	80	100	150',
            'S3	100	15	85	-	G3	3000	300	400	70	100	125',
            'S3	100	20	80	+	G1	1000	400	500	60	100	100',
            'S4	100	25	75	-	G2	2000	500	600	50	100	75',
            'S4	100	30	70	+	G3	3000	600	700	40	100	50',
            '# end'
        )
        obs = list(parse_paf_file_ft(iter(paf), set()))
        exp = [
            ('S1', {'G1'}),
            ('S2', {'G1', 'G2'}),
            ('S3', {'G3', 'G1'}),
            ('S4', {'G2', 'G3'})
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        obs = list(parse_paf_file_ft(iter(paf), {'G1'}))
        self.assertEqual(len(obs), 1)
        self.assertTupleEqual(obs[0], ('S4', {'G2', 'G3'}))
        self.assertFalse(list(parse_paf_file_ft(iter(('hello')), {'G1'})))
        self.assertFalse(list(parse_paf_file_ft(iter(()), {'G1'})))

    def test_parse_paf_file_ex_ft(self):
        paf = (
            'S1	100	0	100	+	G1	1000	0	100	100	100	200',
            'S2	100	5	95	-	G1	1000	100	200	90	95	175',
            'S2	100	10	90	+	G2	2000	200	300	80	90	150',
            'S3	100	15	85	-	G3	3000	300	400	70	85	125',
            'S3	100	20	80	+	G1	1000	400	500	60	80	100',
            'S4	100	25	75	-	G2	2000	500	600	50	75	75',
            'S4	100	30	70	+	G3	3000	600	700	40	70	50',
            '# end'
        )
        obs = list(parse_paf_file_ex_ft(iter(paf), set()))
        exp = [
            ('S1', [('G1', 200, 100,   1, 100)]),
            ('S2', [('G1', 175,  95, 101, 200),
                    ('G2', 150,  90, 201, 300)]),
            ('S3', [('G3', 125,  85, 301, 400),
                    ('G1', 100,  80, 401, 500)]),
            ('S4', [('G2',  75,  75, 501, 600),
                    ('G3',  50,  70, 601, 700)])
        ]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTupleEqual(o, e)
        obs = list(parse_paf_file_ex_ft(iter(paf), {'G1'}))
        self.assertEqual(len(obs), 1)
        self.assertTupleEqual(obs[0], exp[3])
        self.assertFalse(list(parse_paf_file_ex_ft(iter(('hello')), {'G1'})))
        self.assertFalse(list(parse_paf_file_ex_ft(iter(()), {'G1'})))


if __name__ == '__main__':
    main()
