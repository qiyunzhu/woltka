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
from os import listdir
from os.path import isdir, join, basename
from itertools import chain
import click

from .table import (
    read_table, table_shape, table_max_f, frac_table, divide_table,
    scale_table, round_table, write_table, filter_table, merge_tables,
    add_metacol, clip_table, collapse_table, calc_coverage)
from .file import readzip, read_map_1st, read_map_many, read_map_all
from .util import scale_factor
from .tree import read_names
from .ordinal import load_gene_lens


def normalize_wf(input_fp:  str,
                 output_fp: str,
                 sizes_fp:  str = None,
                 scale:     str = None,
                 digits:    int = None):
    """Workflow for normalizing a profile.

    Raises
    ------
    SystemExit
        Scale factor is not valid.
    SystemExit
        Feature not found in size map.
    """
    # read input profile
    table, _ = read_table(input_fp)

    # estimate maximum decimal precision
    if digits is None:
        digits = table_max_f(table)

    # normalize to fractions
    if not sizes_fp:
        click.echo('Normalizing profile to fractions...', nl=False)
        table = frac_table(table)
        click.echo(' Done.')

    # filter profile by threshold
    else:
        click.echo(f'Reading feature sizes from file: {basename(sizes_fp)}...',
                   nl=False)
        with readzip(sizes_fp, {}) as f:

            # check file format based on first line
            try:
                line = next(f)
            except StopIteration:
                exit('Size map file is empty or unreadable.')
            f_ = chain(iter([line]), f)

            # gene coordinates file
            if line[0] in '>#':
                sizes = load_gene_lens(f_)

            # regular mapping file
            else:
                sizes = {k: float(v) for k, v in read_map_1st(f_)}
        click.echo(' Done.')
        click.echo('Normalizing profile by feature size...', nl=False)
        try:
            divide_table(table, sizes)
        except KeyError:
            exit('One or more features are not found in the size map.')
        click.echo(' Done.')

    # scale table values
    if scale:
        try:
            scale = scale_factor(scale)
        except ValueError:
            exit(f'"{scale}" is not a valid scale factor.')
        click.echo(f'Scaling profile by {scale} times...', nl=False)
        scale_table(table, scale)
        click.echo(' Done.')

    # round table values
    round_table(table, digits or None)

    # write normalized profile
    write_table(table, output_fp)
    click.echo('Normalized profile written.')


def filter_wf(input_fp:      str,
              output_fp:     str,
              min_count:     int = None,
              min_percent: float = None):
    """Workflow for filtering a profile based on a per-sample threshold.

    Raises
    ------
    SystemExit
        Neither threshold is provided.
    SystemExit
        Both thresholds are provided.
    SystemExit
        Percentage threshold >= 100.

    See Also
    --------
    .cli.filter_cmd
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

    # read input profile
    table, _ = read_table(input_fp)
    n = table_shape(table)[0]
    click.echo(f'Number of features before filtering: {n}.')

    # filter profile by threshold
    click.echo('Filtered profile...', nl=False)
    table = filter_table(table, th)
    click.echo(' Done.')
    n = table_shape(table)[0]
    click.echo(f'Number of features after filtering: {n}.')

    # write filtered profile
    write_table(table, output_fp)
    click.echo('Filtered profile written.')


def merge_wf(input_fps: list,
             output_fp:  str):
    """Workflow for merging multiple profiles into one.

    Raises
    ------
    SystemExit
        Only one profile is provided.
    SystemExit
        Found invalid profile.

    See Also
    --------
    .cli.merge_cmd
        Command-line arguments and help information.
    """
    click.echo('Reading profiles...')
    tables = []

    def _read_profile(fp):
        try:
            table, _ = read_table(fp)
        except ValueError:
            exit(f'Cannot parse {basename(fp)} as a profile.')
        n, m = table_shape(table)
        click.echo(f'  Read {basename(fp)}. Samples: {m}, features: {n}.')
        tables.append(table)

    for fp in input_fps:
        if isdir(fp):
            for fname in listdir(fp):
                _read_profile(join(fp, fname))
        else:
            _read_profile(fp)

    if len(tables) == 1:
        exit('Please provide two or more profiles.')
    click.echo(f'Done. Number of profiles read: {len(tables)}.')

    # estimate maximum decimal precision
    digits = max(map(table_max_f, tables))

    # merge profiles
    click.echo('Merging profiles...', nl=False)
    table = merge_tables(tables)
    click.echo(' Done.')
    n, m = table_shape(table)
    click.echo(f'Number of samples after merging: {m}.')
    click.echo(f'Number of features after merging: {n}.')

    # round table values
    round_table(table, digits or None)

    # write merged profile
    write_table(table, output_fp)
    click.echo('Merged profile written.')


def collapse_wf(input_fp:  str,
                output_fp: str,
                map_fp:    str = None,
                divide:   bool = False,
                field:     int = None,
                nested:   bool = False,
                sep:       str = None,
                names_fp:  str = None):
    """Workflow for collapsing a profile based on many-to-many mapping.

    Raises
    ------
    SystemExit
        No mapping is found in mapping file.

    See Also
    --------
    .cli.collapse_cmd
        Command-line arguments and help information.
    """
    # read input profile
    table, _ = read_table(input_fp)
    n = table_shape(table)[0]
    click.echo(f'Number of features before collapsing: {n}.')

    # read mapping file (many-to-many okay)
    if map_fp:
        fname = basename(map_fp)
        click.echo(f'Reading mapping file: {fname}...', nl=False)
        with readzip(map_fp, {}) as f:
            mapping = read_map_many(f)
        click.echo(' Done.')
        if not mapping:
            exit(f'No source-target mapping is found in {fname}.')

    # determine default field separator
    if sep is None:
        sep = '_' if nested else '|'

    click.echo('Collapsing profile...', nl=False)

    # maximum decimal precision
    digits = table_max_f(table)

    # collapse profile by mapping
    if map_fp:
        table = collapse_table(table, mapping, divide, field, sep, nested)

    # just clip feature names
    else:
        table = clip_table(table, field, sep, nested)

    # append feature names (optional)
    if names_fp:
        with readzip(names_fp, {}) as f:
            namedic = read_names(f)
        add_metacol(table, namedic, 'Name')

    # round table values
    round_table(table, digits or None)

    click.echo(' Done.')

    n = table_shape(table)[0]
    click.echo(f'Number of features after collapsing: {n}.')

    # write collapsed profile
    write_table(table, output_fp)
    click.echo('Collapsed profile written.')


def coverage_wf(input_fp:   str,
                map_fp:     str,
                output_fp:  str,
                threshold:  int = None,
                count:     bool = False,
                names_fp:   str = None):
    """Calculate coverage of table over feature groups given membership
    information.
    """
    # read input profile
    table, _ = read_table(input_fp)
    n = table_shape(table)[0]
    click.echo(f'Number of features in profile: {n}.')

    # read mapping file (many-to-many okay)
    with readzip(map_fp, {}) as f:
        mapping = dict(read_map_all(f))
    if not mapping:
        exit(f'No group membership is found in {basename(map_fp)}.')

    # calculate group coverage
    click.echo('Calculating coverage...', nl=False)
    table = calc_coverage(table, mapping, threshold, count)
    click.echo(' Done.')
    n = table_shape(table)[0]
    click.echo(f'Number of feature groups with coverage: {n}.')

    # append feature group names (optional)
    if names_fp:
        with readzip(names_fp, {}) as f:
            namedic = read_names(f)
        add_metacol(table, namedic, 'Name')

    # write coverage table
    write_table(table, output_fp)
    click.echo('Coverage table written.')
