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

from woltka.util import readzip
from woltka.align import (
    parse_align_file, Plain, infer_align_format, assign_parser, parse_map_line,
    parse_b6o_line, parse_sam_line, cigar_to_lens)


class AlignTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_parse_align_file(self):
        mapper = Plain()

        # simple case
        aln = ('R1	G1',
               'R2	G1',
               'R2	G2',
               'R3	G1',
               'R3	G3',
               'R4	G4',
               'R5	G5')
        obs = list(parse_align_file(iter(aln), mapper))
        exp = [{'R1': {'G1'},
                'R2': {'G1', 'G2'},
                'R3': {'G1', 'G3'},
                'R4': {'G4'},
                'R5': {'G5'}}]
        self.assertListEqual(obs, exp)

        # chunk of 1
        obs = list(parse_align_file(iter(aln), mapper, n=1))
        exp = [{'R1': {'G1'}},
               {'R2': {'G1', 'G2'}},
               {'R3': {'G1', 'G3'}},
               {'R4': {'G4'}},
               {'R5': {'G5'}}]
        self.assertListEqual(obs, exp)

        # chunk of 2
        obs = list(parse_align_file(iter(aln), mapper, n=2))
        exp = [{'R1': {'G1'}},
               {'R2': {'G1', 'G2'}},
               {'R3': {'G1', 'G3'}},
               {'R4': {'G4'},
                'R5': {'G5'}}]
        self.assertListEqual(obs, exp)

        # chunk of 3
        obs = list(parse_align_file(iter(aln), mapper, n=3))
        exp = [{'R1': {'G1'},
                'R2': {'G1', 'G2'}},
               {'R3': {'G1', 'G3'},
                'R4': {'G4'}},
               {'R5': {'G5'}}]
        self.assertListEqual(obs, exp)

        # chunk of 4
        obs = list(parse_align_file(iter(aln), mapper, n=4))
        exp = [{'R1': {'G1'},
                'R2': {'G1', 'G2'}},
               {'R3': {'G1', 'G3'},
                'R4': {'G4'},
                'R5': {'G5'}}]
        self.assertListEqual(obs, exp)

        # chunk of 5
        obs = list(parse_align_file(iter(aln), mapper, n=5))
        exp = [{'R1': {'G1'},
                'R2': {'G1', 'G2'},
                'R3': {'G1', 'G3'}},
               {'R4': {'G4'},
                'R5': {'G5'}}]
        self.assertListEqual(obs, exp)

        # format is given
        obs = list(parse_align_file(iter(aln), mapper, fmt='map', n=5))
        self.assertListEqual(obs, exp)

        # empty alignment
        obs = list(parse_align_file(iter([]), mapper))
        self.assertListEqual(obs, [])

        # bad alignment
        with self.assertRaises(ValueError) as ctx:
            list(parse_align_file(iter(('Hi there!',)), mapper))
        self.assertEqual(str(ctx.exception), (
            'Cannot determine alignment file format.'))

    def test_plain_parse(self):
        mapper = Plain()
        parser = assign_parser('map')
        self.assertEqual(mapper.parse('R1	G1', parser), 'R1')
        with self.assertRaises(TypeError):
            mapper.parse('Hi there!', parser)

    def test_plain_append(self):
        mapper = Plain()
        parser = assign_parser('map')
        mapper.parse('R1	G1', parser)
        mapper.append()
        self.assertDictEqual(mapper.map, {'R1': ['G1']})
        mapper.parse('R2	G1', parser)
        mapper.append()
        self.assertDictEqual(mapper.map, {'R1': ['G1'], 'R2': ['G1']})
        mapper.parse('R2	G2', parser)
        mapper.append()
        self.assertDictEqual(
            mapper.map, {'R1': ['G1'], 'R2': ['G1', 'G2']})

    def test_plain_flush(self):
        mapper = Plain()
        parser = assign_parser('map')
        aln = ('R1	G1',
               'R2	G1',
               'R2	G2',
               'R3	G1',
               'R3	G3',
               'R4	G4',
               'R5	G5')
        for line in aln:
            mapper.parse(line, parser)
            mapper.append()
        obs = mapper.flush()
        exp = {'R1': {'G1'},
               'R2': {'G1', 'G2'},
               'R3': {'G1', 'G3'},
               'R4': {'G4'},
               'R5': {'G5'}}
        self.assertDictEqual(obs, exp)

    def test_infer_align_format(self):
        # simple cases
        # map
        line = 'S1/1	NC_123456'
        self.assertEqual(infer_align_format(line), 'map')

        # b6o
        line = 'S1/1	NC_123456	100	100	0	0	1	100	25	124	0.1	100'
        self.assertEqual(infer_align_format(line), 'b6o')

        # sam
        line = 'S1	77	NC_123456	26	0	100M	*	0	0	*	*'
        self.assertEqual(infer_align_format(line), 'sam')

        # sam header
        line = '@HD	VN:1.0	SO:unsorted'
        self.assertEqual(infer_align_format(line), 'sam')

        # cannot determine
        line = 'Hi there!'
        with self.assertRaises(ValueError) as ctx:
            infer_align_format(line)
        self.assertEqual(str(ctx.exception), (
            'Cannot determine alignment file format.'))

        # real files
        # Bowtie2 (sam)
        with readzip(join(self.datdir, 'align', 'bowtie2', 'S01.sam.xz')) as f:
            self.assertEqual(infer_align_format(next(f)), 'sam')

        # BURST (b6o)
        with readzip(join(self.datdir, 'align', 'burst', 'S01.b6.bz2')) as f:
            self.assertEqual(infer_align_format(next(f)), 'b6o')

    def test_assign_parser(self):
        self.assertEqual(assign_parser('map'), parse_map_line)
        self.assertEqual(assign_parser('sam'), parse_sam_line)
        self.assertEqual(assign_parser('b6o'), parse_b6o_line)

    def test_parse_map_line(self):
        aln = ('R1	G1', 'R2	G2	# note', '# end of file')

        # first line
        self.assertTupleEqual(parse_map_line(aln[0]), ('R1', 'G1'))

        # second line (with additional column)
        self.assertTupleEqual(parse_map_line(aln[1]), ('R2', 'G2'))

        # third line (not tab-delimited)
        self.assertIsNone(parse_map_line(aln[2]))

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


if __name__ == '__main__':
    main()
