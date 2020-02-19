#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""QIIME 2 plugin for Woltka.

Notes
-----
Minimum requirement: QIIME 2 2019.1
    This is because Woltka is dependent on Python 3.6+, which was adopted by
    QIIME 2 since 2019.1.
"""

import numpy as np
import biom

from woltka.classify import classify, prep_table


def gotu(align_dir: str) -> biom.Table:
    """Generate a gOTU table based on sequence alignments.
    """
    profile = classify(align_dir, None)['none']
    return make_biom(profile)


def filter_values(table: biom.Table, th: float) -> biom.Table:
    """Filter out low-abundance features within each sample in a table.
    """
    def filter_otus(data, id_, md):
        bound = th if th > 1 else data.sum() * th
        data[data < bound] = 0
        return data

    table.transform(filter_otus, axis='sample')
    table.remove_empty(axis='observation')
    return table


def make_biom(profile: dict) -> biom.Table:
    """Convert a profile into a BIOM table.
    """
    data, index, columns = prep_table(profile)
    return biom.Table(np.array(data), index, columns)
