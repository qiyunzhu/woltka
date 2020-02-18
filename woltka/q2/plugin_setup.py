#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, Str
from q2_types.feature_table import FeatureTable, Frequency

from woltka import __version__
from woltka.q2.plugin import gotu


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
        'align_dir': Str,
    },
    parameter_descriptions={
        'align_dir': 'Directory containing input sequence alignments.',
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
