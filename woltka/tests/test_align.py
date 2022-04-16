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
    plain_mapper, range_mapper, infer_align_format, assign_parser,
    parse_map_file, parse_b6o_file, parse_b6o_file_ext, parse_sam_file,
    parse_sam_file_ext, cigar_to_lens, cigar_to_lens_ord, parse_kraken,
    parse_centrifuge, parse_sam_file_pd)


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
               (['R3'], [{'G1', 'G3'}]),
               (['R4', 'R5'], [{'G4'}, {'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 3
        aln.seek(0)
        obs = plain_mapper(aln, n=3)
        exp = ((['R1', 'R2'], [{'G1'}, {'G1', 'G2'}]),
               (['R3', 'R4'], [{'G1', 'G3'}, {'G4'}]),
               (['R5'], [{'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 4
        aln.seek(0)
        obs = plain_mapper(aln, n=4)
        exp = ((['R1', 'R2', 'R3'], [{'G1'}, {'G1', 'G2'}, {'G1', 'G3'}]),
               (['R4', 'R5'], [{'G4'}, {'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # chunk of 5
        aln.seek(0)
        obs = plain_mapper(aln, n=5)
        exp = ((['R1', 'R2', 'R3'], [{'G1'}, {'G1', 'G2'}, {'G1', 'G3'}]),
               (['R4', 'R5'], [{'G4'}, {'G5'}]))
        self.assertTupleEqual(_res2lst(obs), exp)

        # format is given
        aln.seek(0)
        obs = plain_mapper(aln, fmt='map', n=5)
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
            'R4	G2	95	20	0	0	20	1	41	22	1	1',
            'R5	G3	95	20	0	0	20	1	30	49	1	1',
            'R5	G3	95	20	0	0	20	1	50	69	1	1',
            'Rx	Gx	95	20	0	0	1	20	0	0	1	1',
            '# this is not an alignment')))
        obs = _res2lst(range_mapper(aln))[0]
        exp = [('R1', {'G1': [10, 29]}),
               ('R2', {'G1': [16, 35]}),
               ('R3', {'G2': [21, 39]}),
               ('R4', {'G2': [22, 41]}),
               ('R5', {'G3': [30, 49, 50, 69]})]
        self.assertListEqual(list(obs[0]), [x[0] for x in exp])
        self.assertListEqual(list(obs[1]), [x[1] for x in exp])

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
        self.assertEqual(assign_parser('b6o'), parse_b6o_file)
        self.assertEqual(assign_parser('map', True), parse_map_file)
        self.assertEqual(assign_parser('sam', True), parse_sam_file_ext)
        self.assertEqual(assign_parser('b6o', True), parse_b6o_file_ext)
        with self.assertRaises(ValueError) as ctx:
            assign_parser('xyz')
        self.assertEqual(str(ctx.exception), (
            'Invalid format code: "xyz".'))

    def test_parse_map_file(self):
        aln = (
            'R1	G1',          # 1st line (normal)
            'R2	G2	# note',  # 2nd line (with additional column)
            '# end of file')  # 3rd line (not tab-delimited)
        obs = list(parse_map_file(aln))
        self.assertEqual(len(obs), 2)
        self.assertTupleEqual(obs[0], ('R1', 'G1'))
        self.assertTupleEqual(obs[1], ('R2', 'G2'))

    def test_parse_b6o_file(self):
        b6o = (
            '# BLAST result:',
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270')
        obs = list(parse_b6o_file(b6o))
        self.assertEqual(len(obs), 2)
        exp = [('S1/1', 'NC_123456'),
               ('S1/2', 'NC_123456')]
        self.assertTupleEqual(obs[0], exp[0])
        self.assertTupleEqual(obs[1], exp[1])

    def test_parse_b6o_file_ext(self):
        b6o = (
            '# BLAST result:',
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270')
        obs = list(parse_b6o_file_ext(b6o))
        self.assertEqual(len(obs), 2)
        exp = [('S1/1', 'NC_123456', 345, 100, 225, 324),
               ('S1/2', 'NC_123456', 270, 98, 608, 708)]
        self.assertTupleEqual(obs[0], exp[0])
        self.assertTupleEqual(obs[1], exp[1])

    def test_parse_sam_file(self):
        sam = iter((
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            'S2	16	*	0	0	*	*	0	0	*	*'))
        obs = list(parse_sam_file(sam))
        self.assertEqual(len(obs), 3)
        exp = [('S1/1', 'NC_123456'),
               ('S1/2', 'NC_123456'),
               ('S2',   'NC_789012')]
        # 1st line: header
        # 2nd line: normal, fully-aligned, forward strand
        self.assertTupleEqual(obs[0], exp[0])
        # 3rd line: shortened, reverse strand
        self.assertTupleEqual(obs[1], exp[1])
        # 4th line: not perfectly aligned, unpaired
        self.assertTupleEqual(obs[2], exp[2])
        # 5th line: not aligned

        # header only
        sam = iter((
            '@HD	VN:1.0	SO:unsorted',
            '@SQ	SN:G000005825	LN:4249288',
            '@SQ	SN:G000006175	LN:1936387'))
        obs = list(parse_sam_file(sam))
        self.assertEqual(len(obs), 0)

        # file is empty
        obs = list(parse_sam_file(iter(())))
        self.assertEqual(len(obs), 0)

    def test_parse_sam_file_ext(self):
        sam = iter((
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            'S2	16	*	0	0	*	*	0	0	*	*'))
        obs = list(parse_sam_file_ext(sam))
        self.assertEqual(len(obs), 3)
        exp = [('S1/1', 'NC_123456', None, 100, 26, 125),
               ('S1/2', 'NC_123456', None, 80, 151, 230),
               ('S2', 'NC_789012', None, 90, 186, 280)]
        self.assertTupleEqual(obs[0], exp[0])
        self.assertTupleEqual(obs[1], exp[1])
        self.assertTupleEqual(obs[2], exp[2])

        obs = list(parse_sam_file_ext(iter(())))
        self.assertEqual(len(obs), 0)

    def test_cigar_to_lens(self):
        self.assertTupleEqual(cigar_to_lens('150M'), (150, 150))
        self.assertTupleEqual(cigar_to_lens('3M1I3M1D5M'), (11, 12))

    def test_cigar_to_lens_ord(self):
        self.assertTupleEqual(cigar_to_lens_ord('150M'), (150, 150))
        self.assertTupleEqual(cigar_to_lens_ord('3M1I3M1D5M'), (11, 12))

    def test_parse_kraken(self):
        kra = ('C	S1	561	150	561:100 A:10 562:40',
               'U	S2	0	150	1:80 A:40 0:20 A:10')
        obs = list(parse_kraken(kra))
        self.assertEqual(len(obs), 1)
        exp = [('S1', '561')]
        self.assertTupleEqual(obs[0], exp[0])

    def test_parse_centrifuge(self):
        cen = ('readID	seqID	taxID	score	2ndBestScore	'
               'hitLength	queryLength	numMatches',
               'S1	NC_123456	561	125	0	50	150	1')
        obs = list(parse_centrifuge(cen))
        self.assertEqual(len(obs), 1)
        exp = [('S1', 'NC_123456', 125, 50)]
        self.assertTupleEqual(obs[0], exp[0])

    def test_parse_sam_file_pd(self):
        self.assertIsNone(parse_sam_file_pd([]))


if __name__ == '__main__':
    main()
