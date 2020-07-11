#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, Str, MetadataColumn, Categorical

from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy

from woltka import __version__
from woltka.q2.plugin import gotu, classify


plugin = Plugin(
    name='woltka',
    version=__version__,
    website="https://github.com/qiyunzhu/woltka",
    citations=[],
    short_description='Plugin for metagenomics analysis',
    description=('This is a QIIME 2 plugin supporting various analyses of '
                 'shotgun metagenome datasets'),
    package='woltka')


plugin.methods.register_function(
    function=gotu,
    inputs={},
    input_descriptions={},
    parameters={
        'input_path': Str,
    },
    parameter_descriptions={
        'input_path': ('Path to a multiplexed alignment file, or a directory '
                       'of per-sample alignment files.')
    },
    outputs=[
        ('table', FeatureTable[Frequency])
    ],
    output_descriptions={
        'table': 'Output gOTU table.'
    },
    name='gOTU table generation',
    description=('Generate a gOTU table based on sequence alignments against '
                 'a reference genome database.'),
    citations=[]
)


plugin.methods.register_function(
    function=classify,
    inputs={
        'lineage': FeatureData[Taxonomy],
    },
    input_descriptions={
        'lineage': 'Lineage strings. Can accept Greengenes-style rank prefix.'
    },
    parameters={
        'input_path': Str,
        'rank': Str,
        'alignments_mapping': MetadataColumn[Categorical]
    },
    parameter_descriptions={
        'input_path': (
            'Path to a multiplexed alignment file, or a directory of per-'
            'sample alignment files.'),
        'rank': (
            'Classify sequences at this rank. Enter "none" to directly report '
            'subjects; enter "free" for free-rank classification.'),
        'alignments_mapping': 'Sample metadata'
    },
    outputs=[
        ('table', FeatureTable[Frequency])
    ],
    output_descriptions={
        'table': 'Output feature table.'
    },
    name='Flexible hierarchical sequence classifier',
    description=('Classify sequences based on their alignments to references '
                 'through a hierarchical classification system.'),
    citations=[]
)
