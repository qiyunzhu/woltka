#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import main
from os.path import join
from shutil import rmtree
from tempfile import mkdtemp

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError

from woltka.q2 import (
    SeqAlnMapFormat, PaimApFmtFormat, BLAST6OutFormat, SimpleMapFormat,
    NCBINodesFormat, GeneCoordFormat)


class FormatTests(TestPluginBase):
    package = 'woltka.q2.tests'

    def setUp(self):
        self.tmpdir = mkdtemp()

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_seqalnmap_format_validate_positive(self):
        fp = join(self.tmpdir, 'tmp.sam')
        sam = (
            '@HD	VN:1.0	SO:unsorted',
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*',
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*'
        )
        with open(fp, 'w') as fh:
            print(*sam, sep='\n', file=fh)
        SeqAlnMapFormat(fp, mode='r').validate()

    def test_seqalnmap_format_validate_negative(self):
        fp = join(self.tmpdir, 'tmp.sam')
        sam = ('@hello', '@world!')
        with open(fp, 'w') as fh:
            print(*sam, sep='\n', file=fh)
        errmsg = 'No SAM alignment is found in the file.'
        with self.assertRaisesRegex(ValidationError, errmsg):
            SeqAlnMapFormat(fp, mode='r').validate()

        errmsg = 'Invalid SAM alignment:'
        for sam in (
            '	77	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	x	NC_123456	26	0	100M	*	0	0	*	*',
            'S1	77		26	0	100M	*	0	0	*	*',
            'S1	77	NC_123456	x	0	100M	*	0	0	*	*',
            'S1	77	NC_123456	26	0	x	*	0	0	*	*'
        ):
            with open(fp, 'w') as fh:
                print(sam, file=fh)
            with self.assertRaisesRegex(ValidationError, errmsg):
                SeqAlnMapFormat(fp, mode='r').validate()

    def test_paimapfmt_format_validate_positive(self):
        fp = join(self.tmpdir, 'tmp.paf')
        paf = (
            '# Minimap2 result:',
            'S1/1	150	25	100	+	NC_123456	1000000	1025	1100	70	75	20'
            '	tp:A:P	cm:i:10	s1:i:100	s2:i:0	rl:i:0',
            'S1/2	150	25	125	-	NC_123456	1000000	1200	1300	90	100	6'
            '	tp:A:P	cm:i:15	s1:i:150	s2:i:120	rl:i:0',
            'S2/1	250	50	200	-	NC_789012	500000	100000	100150	100	150	65'
            '	tp:A:P	cm:i:24	s1:i:214	s2:i:193	rl:i:0',
            'S2/2	250	0	250	+	NC_789012	500000	99900	100150	245	250	0'
            '	tp:A:P	cm:i:24	s1:i:282	s2:i:267	rl:i:0'
        )
        with open(fp, 'w') as fh:
            print(*paf, sep='\n', file=fh)
        PaimApFmtFormat(fp, mode='r').validate()

    def test_paimapfmt_format_validate_negative(self):
        fp = join(self.tmpdir, 'tmp.paf')
        paf = ('#hello', '#world!')
        with open(fp, 'w') as fh:
            print(*paf, sep='\n', file=fh)
        errmsg = 'No PAF alignment is found in the file.'
        with self.assertRaisesRegex(ValidationError, errmsg):
            PaimApFmtFormat(fp, mode='r').validate()

        errmsg = 'Invalid PAF alignment:'
        for paf in (
            '	150	25	100	+	NC_123456	1000000	1025	1100	70	75	20',
            'S1	150	25	100	+		1000000	1025	1100	70	75	20',
            'S1	150	25	100	+	NC_123456	1000000	x	1100	70	75	20',
            'S1	150	25	100	+	NC_123456	1000000	1025	x	70	75	20',
            'S1	150	25	100	+	NC_123456	1000000	9999	1100	70	75	20',
            'S1	150	25	100	+	NC_123456	1000000	1025	1100	70	0	20',
            'S1	150	25	100	+	NC_123456	1000000	1025	1100	70	75	x'
        ):
            with open(fp, 'w') as fh:
                print(paf, file=fh)
            with self.assertRaisesRegex(ValidationError, errmsg):
                PaimApFmtFormat(fp, mode='r').validate()

    def test_blast6out_format_validate_positive(self):
        fp = join(self.tmpdir, 'tmp.b6o')
        b6o = (
            '# BLAST result:',
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/2	NC_123456	95	98	2	1	2	99	708	608	3.4e-20	270'
        )
        with open(fp, 'w') as fh:
            print(*b6o, sep='\n', file=fh)
        BLAST6OutFormat(fp, mode='r').validate()

    def test_blast6out_format_validate_negative(self):
        fp = join(self.tmpdir, 'tmp.b6o')
        b6o = ('#hello', '#world!')
        with open(fp, 'w') as fh:
            print(*b6o, sep='\n', file=fh)
        errmsg = 'No BLAST hit record is found in the file.'
        with self.assertRaisesRegex(ValidationError, errmsg):
            BLAST6OutFormat(fp, mode='r').validate()

        errmsg = 'Invalid BLAST hit record:'
        for b6o in (
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30',
            '	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/1		100	100	0	0	1	100	225	324	1.2e-30	345',
            'S1/1	NC_123456	100	0	0	0	1	100	225	324	1.2e-30	345',
            'S1/1	NC_123456	100	100	0	0	1	100	x	324	1.2e-30	345',
            'S1/1	NC_123456	100	100	0	0	1	100	225	x	1.2e-30	345',
            'S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	x'
        ):
            with open(fp, 'w') as fh:
                print(b6o, file=fh)
            with self.assertRaisesRegex(ValidationError, errmsg):
                BLAST6OutFormat(fp, mode='r').validate()

    def test_simplemap_format_validate_positive(self):
        fp = join(self.tmpdir, 'tmp.aln')
        aln = (
            'R1	G1',
            'R2	G2	# note',
            '# end of file'
        )
        with open(fp, 'w') as fh:
            print(*aln, sep='\n', file=fh)
        SimpleMapFormat(fp, mode='r').validate()

    def test_simplemap_format_validate_negative(self):
        fp = join(self.tmpdir, 'tmp.aln')
        aln = (
            '#start',
            'hello',
            'world!',
            'R1		G1',
            '	G1',
            'R1	',
            '-end-'
        )
        with open(fp, 'w') as fh:
            print(*aln, sep='\n', file=fh)
        errmsg = 'No mapping is found in the file.'
        with self.assertRaisesRegex(ValidationError, errmsg):
            SimpleMapFormat(fp, mode='r').validate()

    def test_ncbinodes_format_validate_positive(self):
        fp = join(self.tmpdir, 'tmp.dmp')
        dmp = ('1	|	1	|	no rank	|',
               '2	|	1	|	superkingdom	|',
               '1224	|	2	|	phylum	|',
               '1236	|	1224	|	class	|',
               '32066	|	2	|	phylum	|')
        with open(fp, 'w') as fh:
            print(*dmp, sep='\n', file=fh)
        NCBINodesFormat(fp, mode='r').validate()

        dmp = ('1	1	no rank',
               '2	1	superkingdom',
               '1224	2	phylum',
               '1236	1224	class',
               '32066	2	phylum')
        with open(fp, 'w') as fh:
            print(*dmp, sep='\n', file=fh)
        NCBINodesFormat(fp, mode='r').validate()

    def test_ncbinodes_format_validate_negative(self):
        fp = join(self.tmpdir, 'tmp.dmp')
        errmsg = 'Invalid node mapping:'
        for dmp in ('hello', '	2', '1	', '1|2'):
            with open(fp, 'w') as fh:
                print(dmp, file=fh)
            with self.assertRaisesRegex(ValidationError, errmsg):
                NCBINodesFormat(fp, mode='r').validate()

    def test_genecoord_format_validate_positive(self):
        fp = join(self.tmpdir, 'tmp.coords')
        coords = (
            '>n1',
            'g1	5	29',
            'g2	33	61',
            'g3	65	94'
        )
        with open(fp, 'w') as fh:
            print(*coords, sep='\n', file=fh)
        GeneCoordFormat(fp, mode='r').validate()

        coords = (
            '## GCF_000123456',
            '# NC_123456',
            '1	5	384',
            '2	410	933',
            '# NC_789012',
            '1	912	638',
            '2	529	75'
        )
        with open(fp, 'w') as fh:
            print(*coords, sep='\n', file=fh)
        GeneCoordFormat(fp, mode='r').validate()

    def test_genecoord_format_validate_negative(self):
        fp = join(self.tmpdir, 'tmp.coords')

        errmsg = 'No gene coordinates are found in the file.'
        for coord in ('', '#hello', 'hello	5	10'):
            with open(fp, 'w') as fh:
                print(coord, file=fh)
            with self.assertRaisesRegex(ValidationError, errmsg):
                GeneCoordFormat(fp, mode='r').validate()

        errmsg = 'Cannot extract coordinates from line:'
        for coord in ('hi	x	10', 'hi	5	x', 'hi	5', 'hi'):
            with open(fp, 'w') as fh:
                print(*('#hello', coord), sep='\n', file=fh)
            with self.assertRaisesRegex(ValidationError, errmsg):
                GeneCoordFormat(fp, mode='r').validate()


if __name__ == '__main__':
    main()
