#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import main
from qiime2.plugin.testing import TestPluginBase

from woltka.q2 import (
    SeqAlnMap, PaimApFmt, BLAST6Out, SimpleMap, NCBINodes, GeneCoordinates)


class TypeTests(TestPluginBase):
    package = 'woltka.q2.tests'

    def test_seqalnmap_semantic_type_registration(self):
        self.assertRegisteredSemanticType(SeqAlnMap)

    def test_paimapfmt_semantic_type_registration(self):
        self.assertRegisteredSemanticType(PaimApFmt)

    def test_blast6out_semantic_type_registration(self):
        self.assertRegisteredSemanticType(BLAST6Out)

    def test_simplemap_semantic_type_registration(self):
        self.assertRegisteredSemanticType(SimpleMap)

    def test_ncbinodes_semantic_type_registration(self):
        self.assertRegisteredSemanticType(NCBINodes)

    def test_genecoordinates_semantic_type_registration(self):
        self.assertRegisteredSemanticType(GeneCoordinates)


if __name__ == '__main__':
    main()
