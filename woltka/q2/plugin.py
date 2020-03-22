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

import biom

from woltka.biom import profile_to_biom
from woltka.workflow import workflow


def gotu(input_path: str) -> biom.Table:
    """Generate a gOTU table based on sequence alignments.
    """
    profile = workflow(input_path, None)['none']
    return profile_to_biom(profile)


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
