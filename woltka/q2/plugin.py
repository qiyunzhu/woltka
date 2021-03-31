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

import biom
from pandas import Series
from skbio import TreeNode

from woltka.workflow import build_mapper, classify as cwf
from woltka.file import read_map_1st, read_map_many, read_map_all
from woltka.table import prep_table, calc_coverage
from woltka.biom import table_to_biom, filter_biom, collapse_biom
from woltka.tree import read_nodes, read_lineage, read_newick, fill_root
from woltka.__init__ import __name__, __version__


def classify(alignment:             str,
             target_rank:           str,
             reference_taxonomy: Series = None,
             reference_tree:   TreeNode = None,
             reference_nodes:       str = None,
             taxon_map:             str = None,
             trim_subject:         bool = False,
             gene_coordinates:      str = None,
             overlap_threshold:     int = 80,
             unique_assignment:    bool = False,
             majority_threshold:    int = None,
             above_given_rank:     bool = False,
             subject_is_okay:      bool = False,
             report_unassigned:    bool = False) -> biom.Table:
    """Classify sequences based on their alignments to references through a
    hierarchical classification system.
    """
    # validate classification system
    num_ref = len(list(filter(None.__ne__, (
        reference_taxonomy, reference_tree, reference_nodes))))
    if num_ref > 1:
        raise ValueError('Only one reference classification system can be '
                         'specified.')
    elif num_ref == 0 and target_rank != 'none':
        raise ValueError('A reference classification system must be specified '
                         f'for classification at the rank "{target_rank}".')

    # build classification hierarchy
    tree, rankdic, namedic = {}, {}, {}

    # read taxonomy
    if reference_taxonomy is not None:
        tree, rankdic = read_lineage(StringIO(reference_taxonomy.to_csv(
            sep='\t', header=False)))

    # read phylogeny
    if reference_tree is not None:
        tree = read_newick(StringIO(str(reference_tree)))

    # read taxdump
    if reference_nodes is not None:
        with open(reference_nodes, 'r') as fh:
            tree, rankdic = read_nodes(fh)

    # read taxon mapping
    if taxon_map is not None:
        with open(taxon_map, 'r') as fh:
            tree.update(read_map_1st(fh))

    # fill root
    root = fill_root(tree)

    # build mapping module
    mapper, chunk = build_mapper(
        coords_fp=gene_coordinates, overlap=overlap_threshold)

    # classify query sequences
    profile = cwf(mapper=mapper,
                  files=[alignment],
                  demux=True,
                  trimsub=trim_subject and '_',
                  tree=tree,
                  rankdic=rankdic,
                  namedic=namedic,
                  root=root,
                  ranks=[target_rank],
                  uniq=unique_assignment,
                  major=majority_threshold,
                  above=above_given_rank,
                  subok=subject_is_okay,
                  unasgd=report_unassigned,
                  chunk=chunk,
                  zippers={})[target_rank]

    # generate feature table
    table = table_to_biom(*prep_table(
        profile, rankdic=rankdic, namedic=namedic))
    table.generated_by = f'{__name__}-{__version__}'

    return table


def psfilter(table:  biom.Table,
             min_count:     int = None,
             min_percent: float = None) -> biom.Table:
    """Filter a feature table by per-sample abundance of features.
    """
    # validate parameters
    if not any((min_count, min_percent)):
        raise ValueError('Please specify either minimum count or minimum '
                         'percentage threshold.')
    if all((min_count, min_percent)):
        raise ValueError('Only one of minimum count or minimum percentage '
                         'thresholds can be specified.')
    if min_percent and min_percent >= 100:
        raise ValueError('Minimum percentage threshold must be below 100.')

    # determine threshold
    th = min_count or min_percent / 100

    # filter feature table
    table = filter_biom(table, th)
    table.generated_by = f'{__name__}-{__version__}'

    return table


def collapse(table: biom.Table,
             mapping:      str,
             normalize:   bool = False) -> biom.Table:
    """Collapse a feature table based on many-to-many mapping.
    """
    with open(mapping, 'r') as fh:
        mapping = read_map_many(fh)
    table = collapse_biom(table, mapping, normalize)
    table.generated_by = f'{__name__}-{__version__}'
    return table


def coverage(table: biom.Table,
             mapping:      str,
             threshold:    int = None,
             count:       bool = False) -> biom.Table:
    """Calculate a feature table's coverage over feature groups.
    """
    with open(mapping, 'r') as fh:
        mapping = dict(read_map_all(fh))
    table = calc_coverage(table, mapping, threshold, count)
    table = table_to_biom(*table)
    table.generated_by = f'{__name__}-{__version__}'
    return table
