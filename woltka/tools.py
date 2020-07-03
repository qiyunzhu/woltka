#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Features under the `tools` command group.
"""

from sys import exit
import click

from .table import read_table, table_shape, filter_table, write_table


def filter_wf(input_fp:      str,
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
        exit('Please specify either minimum count or minimum percentage '
             'threshold.')
    if all((min_count, min_percent)):
        exit('Only one of minimum count or minimum percentage thresholds '
             'can be specified.')
    if min_percent and min_percent >= 100:
        exit('Minimum percentage threshold must be below 100.')

    # determine threshold
    th = min_count or min_percent / 100

    # read input table
    table, fmt = read_table(input_fp)
    n = table_shape(table)[0]
    click.echo(f'Number of features before filtering: {n}.')

    # filter table by threshold
    click.echo(f'Filtered profile...', nl=False)
    table = filter_table(table, th)
    click.echo(' Done.')
    n = table_shape(table)[0]
    click.echo(f'Number of features After filtering: {n}.')

    # write filtered table
    write_table(table, output_fp)
    click.echo('Filtered profile written.', err=True)
