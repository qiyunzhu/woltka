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

from woltka.align import (
    infer_align_format, parse_line_ordinal, parse_b6o_line, parse_sam_line,
    cigar_to_lens, strip_index, demultiplex)


class AlignTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_infer_align_format(self):
        # simple cases
        # map
        line = 'S1/1	NC_123456'
        self.assertEqual(infer_align_format(line), 'map')

        # b6o
        line = 'S1/1	NC_123456	100	100	0	0	1	100	25	124	0.1	100'
        self.assertEqual(infer_align_format(line), 'b6o')

        # sam
        line = '@HD	VN:1.0	SO:unsorted'
        self.assertEqual(infer_align_format(line), 'sam')
        line = 'S1	77	NC_123456	26	0	100M	*	0	0	*	*'
        self.assertEqual(infer_align_format(line), 'sam')

    def test_parse_line_ordinal(self):
        # b6o (BLAST, DIAMOND, BURST, etc.)
        parser = parse_b6o_line
        b6o = ('S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
               'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270')
        rids, lenmap, locmap = [], {}, {}

        # first line
        parse_line_ordinal(b6o[0], parser, rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1'])
        self.assertDictEqual(lenmap, {'NC_123456': {0: 100}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (225, True, False, 0), (324, False, False, 0)]})

        # second line
        parse_line_ordinal(b6o[1], parser, rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1', 'S1/2'])
        self.assertDictEqual(lenmap, {'NC_123456': {0: 100, 1: 98}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (225, True, False, 0), (324, False, False, 0),
            (608, True, False, 1), (708, False, False, 1)]})

        # sam (Bowtie2, BWA, etc.)
        parser = parse_sam_line
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            'S2	16	*	0	0	*	*	0	0	*	*')
        rids, lenmap, locmap = [], {}, {}

        # header
        parse_line_ordinal(sam[0], parser, rids, lenmap, locmap)
        self.assertFalse(len(rids) + len(lenmap) + len(locmap))

        # normal, fully-aligned, forward strand
        parse_line_ordinal(sam[1], parser, rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1'])
        self.assertDictEqual(lenmap, {'NC_123456': {0: 100}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (26, True, False, 0), (125, False, False, 0)]})

        # shortened, reverse strand
        parse_line_ordinal(sam[2], parser, rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1', 'S1/2'])
        self.assertDictEqual(lenmap, {'NC_123456': {
            0: 100, 1: 80}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (26,  True, False, 0), (125, False, False, 0),
            (151, True, False, 1), (230, False, False, 1)]})

        # not perfectly aligned, unpaired
        parse_line_ordinal(sam[3], parser, rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1', 'S1/2', 'S2'])
        self.assertEqual(lenmap['NC_789012'][2], 90)
        self.assertTupleEqual(locmap['NC_789012'][0], (186,  True, False, 2))
        self.assertTupleEqual(locmap['NC_789012'][1], (280, False, False, 2))

        # not aligned
        parse_line_ordinal(sam[4], parser, rids, lenmap, locmap)
        self.assertEqual(len(rids), 3)
        self.assertEqual(len(lenmap['NC_789012']), 1)
        self.assertEqual(len(locmap['NC_789012']), 2)

    def test_parse_b6o_line(self):
        b6o = ('S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
               'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270')

        # first line
        obs = parse_b6o_line(b6o[0])
        exp = ('S1/1', 'NC_123456', 345, 100, 225, 324)
        self.assertTupleEqual(obs, exp)

        # second line
        obs = parse_b6o_line(b6o[1])
        exp = ('S1/2', 'NC_123456', 270, 98, 608, 708)
        self.assertTupleEqual(obs, exp)

    def test_parse_sam_line(self):
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*',
            'S2	16	*	0	0	*	*	0	0	*	*')

        # header
        self.assertIsNone(parse_sam_line(sam[0]))

        # normal, fully-aligned, forward strand
        obs = parse_sam_line(sam[1])
        exp = ('S1/1', 'NC_123456', None, 100, 26, 125)
        self.assertTupleEqual(obs, exp)

        # shortened, reverse strand
        obs = parse_sam_line(sam[2])
        exp = ('S1/2', 'NC_123456', None, 80, 151, 230)
        self.assertTupleEqual(obs, exp)

        # not perfectly aligned, unpaired
        obs = parse_sam_line(sam[3])
        exp = ('S2', 'NC_789012', None, 90, 186, 280)
        self.assertTupleEqual(obs, exp)

        # not aligned
        self.assertIsNone(parse_sam_line(sam[4]))

    def test_cigar_to_lens(self):
        self.assertTupleEqual(cigar_to_lens('150M'), (150, 150))
        self.assertTupleEqual(cigar_to_lens('3M1I3M1D5M'), (11, 12))

        with self.assertRaises(ValueError) as context:
            cigar_to_lens('*')
        msg = 'Missing CIGAR string.'
        self.assertEqual(str(context.exception), msg)

    def test_strip_index(self):
        dic = {'R1': ['G1_1', 'G1_2', 'G2_3', 'G3'],
               'R2': ['G1_1', 'G1.3', 'G4_5', 'G4_x']}
        strip_index(dic)
        self.assertDictEqual(dic, {
            'R1': ['G1', 'G1',   'G2', 'G3'],
            'R2': ['G1', 'G1.3', 'G4', 'G4']})

    def test_demultiplex(self):
        # simple case
        dic = {'S1_R1': 5,
               'S1_R2': 12,
               'S1_R3': 3,
               'S2_R1': 10,
               'S2_R2': 8,
               'S2_R4': 7,
               'S3_R2': 15,
               'S3_R3': 1,
               'S3_R4': 5}
        obs = demultiplex(dic)
        exp = {'S1': {'R1': 5, 'R2': 12, 'R3': 3},
               'S2': {'R1': 10, 'R2': 8, 'R4': 7},
               'S3': {'R2': 15, 'R3': 1, 'R4': 5}}
        self.assertDictEqual(obs, exp)

        # change separator, no result
        obs = demultiplex(dic, sep='.')
        self.assertDictEqual(obs, {'': dic})

        # enforce sample Ids
        obs = demultiplex(dic, samples=['S1', 'S2', 'SX'])
        exp = {x: exp[x] for x in ['S1', 'S2']}
        self.assertDictEqual(obs, exp)


if __name__ == '__main__':
    main()
