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

from io import StringIO

import pandas as pd
import biom
from skbio import TreeNode

from woltka.workflow import (
    workflow, parse_samples, build_mapper, classify as cwf)
from woltka.table import prep_table
from woltka.biom import table_to_biom
from woltka.tree import read_lineage, fill_root


def gotu(input_path: str) -> biom.Table:
    """Generate a gOTU table based on sequence alignments.
    """
    profile = workflow(input_path, None)['none']
    return table_to_biom(*prep_table(profile))


def classify(input_path: str,
             rank: str,
             lineage: pd.Series = None,
             sample_metadata: pd.DataFrame = None) -> biom.Table:
    """Classify sequences based on their alignments to references through a
    hierarchical classification system.
    """
    # available external compressors
    zippers = {}

    # parse input samples
    samples, files, demux = parse_samples(input_path)

    # build classification hierarchy
    tree, rankdic, namedic, root = build_hierarchy(taxonomy=lineage)

    # build mapping module
    mapper, chunk = build_mapper()

    # classify query sequences
    profile = cwf(mapper=mapper,
                  files=files,
                  samples=samples,
                  demux=demux,
                  tree=tree,
                  rankdic=rankdic,
                  namedic=namedic,
                  root=root,
                  ranks=[rank],
                  chunk=chunk,
                  zippers=zippers)[rank]

    return table_to_biom(*prep_table(profile))


def build_hierarchy(taxonomy: pd.Series = None,
                    phylogeny: TreeNode = None,
                    mapping: pd.Series = None) -> (dict, dict, dict, dict):
    tree, rankdic, namedic = {}, {}, {}

    if taxonomy is not None:
        tree, rankdic = read_lineage(StringIO(taxonomy.to_csv(
            sep='\t', header=False)))

    root = fill_root(tree)
    return tree, rankdic, namedic, root


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
