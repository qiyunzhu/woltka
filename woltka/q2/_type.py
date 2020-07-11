#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType
from q2_types.feature_data import FeatureData


SeqAlnMap = SemanticType('SeqAlnMap', variant_of=FeatureData.field['type'])
BLAST6Out = SemanticType('BLAST6Out', variant_of=FeatureData.field['type'])
SimpleMap = SemanticType('SimpleMap', variant_of=FeatureData.field['type'])
NCBINodes = SemanticType('NCBINodes', variant_of=FeatureData.field['type'])
GeneCoordinates = SemanticType('GeneCoordinates')
