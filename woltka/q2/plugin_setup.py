#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from importlib import import_module

from qiime2.plugin import Plugin, Str, Bool, Int, Float, Range
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.tree import Phylogeny, Rooted

from ._format import (SeqAlnMapFormat, SeqAlnMapDirFmt,
                      BLAST6OutFormat, BLAST6OutDirFmt,
                      SimpleMapFormat, SimpleMapDirFmt,
                      NCBINodesFormat, NCBINodesDirFmt,
                      GeneCoordFormat, GeneCoordDirFmt)
from ._type import SeqAlnMap, BLAST6Out, SimpleMap, NCBINodes, GeneCoordinates

from woltka import __version__
from woltka.q2.plugin import classify, psfilter, collapse, coverage


plugin = Plugin(
    name='woltka',
    version=__version__,
    website="https://github.com/qiyunzhu/woltka",
    citations=[],
    short_description='Plugin for metagenomics analysis',
    description=('This is a QIIME 2 plugin supporting various analyses of '
                 'shotgun metagenome datasets'),
    package='woltka')


plugin.register_semantic_types(
    SeqAlnMap, BLAST6Out, SimpleMap, NCBINodes, GeneCoordinates)

plugin.register_formats(SeqAlnMapFormat, SeqAlnMapDirFmt,
                        BLAST6OutFormat, BLAST6OutDirFmt,
                        SimpleMapFormat, SimpleMapDirFmt,
                        NCBINodesFormat, NCBINodesDirFmt,
                        GeneCoordFormat, GeneCoordDirFmt)

plugin.register_semantic_type_to_format(
    FeatureData[SeqAlnMap], artifact_format=SeqAlnMapDirFmt)
plugin.register_semantic_type_to_format(
    FeatureData[BLAST6Out], artifact_format=BLAST6OutDirFmt)
plugin.register_semantic_type_to_format(
    FeatureData[SimpleMap], artifact_format=SimpleMapDirFmt)
plugin.register_semantic_type_to_format(
    FeatureData[NCBINodes], artifact_format=NCBINodesDirFmt)
plugin.register_semantic_type_to_format(
    GeneCoordinates, artifact_format=GeneCoordDirFmt)

alnfmts = SeqAlnMap | BLAST6Out | SimpleMap

plugin.methods.register_function(
    function=classify,
    inputs={
        'alignment': FeatureData[alnfmts],
        'reference_taxonomy': FeatureData[Taxonomy],
        'reference_tree': Phylogeny[Rooted],
        'reference_nodes': FeatureData[NCBINodes],
        'taxon_map': FeatureData[SimpleMap],
        'gene_coordinates': GeneCoordinates
    },
    input_descriptions={
        'alignment': (
            'Multiplexed sequence alignment map to be classified. Can accept '
            'SAM, BLAST6 or simple map format.'),
        'reference_taxonomy': (
            'Reference taxonomic lineage strings.'),
        'reference_tree': (
            'Reference phylogenetic tree.'),
        'reference_nodes': (
            'Reference taxonomic nodes.'),
        'taxon_map': (
            'Mapping of subject IDs to taxon IDs.'),
        'gene_coordinates': (
            'Mapping of gene IDs to start and end positions on host genomes.')
    },
    parameters={
        'target_rank': Str,
        'overlap_threshold': Int % Range(1, 100),
        'trim_subject': Bool,
        'unique_assignment': Bool,
        'majority_threshold': Int % Range(51, 99),
        'above_given_rank': Bool,
        'subject_is_okay': Bool,
        'report_unassigned': Bool
    },
    parameter_descriptions={
        'target_rank': (
            'Classify sequences at this rank. Enter "none" to directly report '
            'subjects; enter "free" for free-rank classification.'),
        'overlap_threshold': (
            'Read/gene overlapping percentage threshold.'),
        'trim_subject': (
            'Trim subject IDs at the last underscore.'),
        'unique_assignment': (
            'One sequence can only be assigned to one classification unit, or '
            'remain unassigned if there is ambiguity. Otherwise, all candidate'
            ' units are reported and their counts are normalized.'),
        'majority_threshold': (
            'In given-rank classification, use majority rule at this '
            'percentage threshold to determine assignment when there are '
            'multiple candidates.'),
        'above_given_rank': (
            'In given-rank classification, allow assigning a sequence to '
            'a higher rank if it cannot be assigned to the current rank.'),
        'subject_is_okay': (
            'In free-rank classification, allow assigning a sequence to its '
            'direct subject, if applicable, before going up in hierarchy.'),
        'report_unassigned': (
            'Report Frequency of unassigned sequences (will be marked as '
            '"Unassigned").')
    },
    outputs=[
        ('classified_table', FeatureTable[Frequency])
    ],
    output_descriptions={
        'classified_table': (
            'The resulting table of frequencies of classification units.')
    },
    name='Flexible hierarchical sequence classifier',
    description=('Classify sequences based on their alignments to references '
                 'through a hierarchical classification system.'),
    citations=[]
)


plugin.methods.register_function(
    function=psfilter,
    inputs={'table': FeatureTable[Frequency]},
    input_descriptions={'table': 'Feature table to be filtered.'},
    parameters={'min_count': Int % Range(1, None),
                'min_percent': Float % Range(0.0, 100)},
    parameter_descriptions={
        'min_count': 'Per-sample minimum count threshold.',
        'min_percent': 'Per-sample minimum count threshold.'},
    outputs=[('filtered_table', FeatureTable[Frequency])],
    output_descriptions={'filtered_table': 'Filtered feature table.'},
    name='Per-sample feature filter',
    description=('Filter a feature table by per-sample feature abundance.'),
    citations=[]
)


plugin.methods.register_function(
    function=collapse,
    inputs={'table': FeatureTable[Frequency],
            'mapping': FeatureData[SimpleMap]},
    input_descriptions={'table': 'Feature table to be collapsed.',
                        'mapping': 'Mapping of source features to targets.'},
    parameters={'normalize': Bool},
    parameter_descriptions={'normalize': (
        'Normalize target counts by number of targets per source.')},
    outputs=[('collapsed_table', FeatureTable[Frequency])],
    output_descriptions={'collapsed_table': 'Collapsed feature table.'},
    name='Many-to-many feature mapper',
    description=('Collapse a feature table based on many-to-many feature '
                 'mapping.'),
    citations=[]
)


plugin.methods.register_function(
    function=coverage,
    inputs={'table': FeatureTable[Frequency],
            'mapping': FeatureData[SimpleMap]},
    input_descriptions={'table': 'Feature table to calculate coverage.',
                        'mapping': 'Mapping of feature groups to members.'},
    parameters={'threshold': Int % Range(1, 100),
                'count': Bool},
    parameter_descriptions={
        'threshold': ('Convert coverage to presence (1) / absence (0) data by '
                      'this percentage threshold.'),
        'count': ('Record numbers of covered features instead of percentages '
                  '(overrides threshold).')},
    outputs=[('coverage_table', FeatureTable[Frequency])],
    output_descriptions={'coverage_table': 'Feature group coverage table.'},
    name='Group coverage calculator',
    description='Calculate per-sample coverage of feature groups.',
    citations=[]
)


import_module('woltka.q2._transformer')
