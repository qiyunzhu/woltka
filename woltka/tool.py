#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Features other than the main `classify` subcommand.
"""

import click
from biom import load_table

from .table import read_tsv, filter_table, write_tsv
from .biom import filter_biom, write_biom, table_to_biom


def filter_table_wf(input_fp:      str,
                    output_fp:     str,
                    min_count:     int = None,
                    min_percent: float = None):
    """Workflow for filtering a feature table based on a per-sample threshold.

    See Also
    --------
    .cli.filter_table
        Command-line arguments and help information.
    """
    # validate parameters
    if not any((min_count, min_percent)):
        click.echo('Please specify either minimum count or minimum percentage '
                   'threshold.', err=True)
    if all((min_count, min_percent)):
        click.echo('Only one of minimum count or minimum percentage thresholds'
                   'can be specified', err=True)
    if min_percent and min_percent >= 100:
        click.echo('Minimum percentage threshold must be below 100.', err=True)

    # determine threshold
    th = min_count or min_percent / 100

    # read input table
    table, fmt = read_table(input_fp)

    # filter table by threshold
    filter_func = filter_table if fmt == 'tsv' else filter_biom
    table = filter_func(table, th)

    # write filtered table
    if output_fp.endswith('.biom'):
        write_biom(table, output_fp)
    else:
        with open(output_fp, 'w') as fh:
            write_tsv(fh, *table)

