#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Standalone command-line interface (CLI) of the program.
"""

import click

from . import __version__
from .workflow import workflow
from .tools import filter_wf, merge_wf, collapse_wf, coverage_wf


class NaturalOrderGroup(click.Group):
    """Natural ordering of click command groups.
    """
    def list_commands(self, ctx):
        return self.commands.keys()


@click.version_option(__version__)
@click.group(cls=NaturalOrderGroup)
def cli():
    pass  # pragma: no cover


# `classify` invokes the main classification workflow

@cli.command('classify')
# input and output
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    help='Path to input alignment file or directory of alignment files.')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True),
    help='Path to output profile file or directory of profile files.')
# input files
@click.option(
    '--format', '-f', 'input_fmt',
    type=click.Choice(['b6o', 'sam', 'map'], case_sensitive=False),
    help=('Format of read alignments: "b6o": BLAST tabular format, "sam": SAM '
          'format, "map": simple map of query <tab> subject. If not specified,'
          ' program will automatically infer from file content.'))
@click.option(
    '--filext', '-e', 'input_ext',
    help='Input filename extension following sample ID.')
@click.option(
    '--samples', '-s', type=click.STRING,
    help='List of sample IDs to be included.')
@click.option(
    '--demux/--no-demux', default=None,
    help='Demultiplex alignment by first underscore in query identifier.')
@click.option(
    '--trim-sub', 'trimsub',
    help='Trim subject IDs at the last given delimiter.')
# hierarchies
@click.option(
    '--nodes', 'nodes_fps', type=click.Path(exists=True), multiple=True,
    help='Hierarchies defined by NCBI nodes.dmp or compatible formats.')
@click.option(
    '--newick', 'newick_fps', type=click.Path(exists=True), multiple=True,
    help='Hierarchies defined by a tree in Newick format.')
@click.option(
    '--lineage', 'lineage_fps', type=click.Path(exists=True), multiple=True,
    help='Lineage strings. Can accept Greengenes-style rank prefix.')
@click.option(
    '--columns', 'columns_fps', type=click.Path(exists=True), multiple=True,
    help='Table of classification units per rank (column).')
@click.option(
    '--map', '-m', 'map_fps', type=click.Path(exists=True), multiple=True,
    help='Mapping of lower classification units to higher ones.')
@click.option(
    '--map-as-rank', is_flag=True,
    help='Extract rank name from map filename.')
@click.option(
    '--names', '-n', 'names_fps', type=click.Path(exists=True), multiple=True,
    help=('Names of classification units as defined by NCBI names.dmp or a '
          'simple map.'))
# assignment
@click.option(
    '--rank', '-r', 'ranks', type=click.STRING,
    help=('Classify sequences at this rank. Enter "none" to directly report '
          'subjects; enter "free" for free-rank classification. Can '
          'specify multiple comma-separated ranks.'))
@click.option(
    '--uniq', is_flag=True,
    help=('One sequence can only be assigned to one classification unit, or '
          'remain unassigned if there is ambiguity. Otherwise, all candidate '
          'units are reported and their counts are normalized.'))
@click.option(
    '--major', type=click.IntRange(51, 99),
    help=('In given-rank classification, use majority rule at this percentage '
          'threshold to determine assignment when there are multiple '
          'candidates.'))
@click.option(
    '--above', is_flag=True,
    help=('In given-rank classification, allow assigning a sequence to '
          'a higher rank if it cannot be assigned to the current rank.'))
@click.option(
    '--subok', is_flag=True,
    help=('In free-rank classification, allow assigning a sequence to its '
          'direct subject, if applicable, before going up in hierarchy.'))
# gene matching
@click.option(
    '--coords', '-c', 'coords_fp', type=click.Path(exists=True),
    help='Reference gene coordinates on genomes.')
@click.option(
    '--overlap', type=click.IntRange(1, 100), default=80, show_default=True,
    help='Read/gene overlapping percentage threshold.')
# stratification
@click.option(
    '--stratify', '-t', 'strata_dir',
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    help='Directory of read-to-feature maps for stratification.')
# output files
@click.option(
    '--to-biom/--to-tsv', 'output_fmt', default=None,
    help='Output profile format (BIOM or TSV).')
@click.option(
    '--unassigned', is_flag=True,
    help='Report unassigned sequences.')
@click.option(
    '--name-as-id', is_flag=True,
    help='Replace feature IDs with names.')
@click.option(
    '--add-rank', is_flag=True,
    help='Append feature ranks to table.')
@click.option(
    '--add-lineage', is_flag=True,
    help='Append lineage strings to table.')
@click.option(
    '--outmap', '-u', 'outmap_dir',
    type=click.Path(dir_okay=True),
    help='Write read-to-feature maps to this directory.')
@click.option(
    '--zipmap', 'outmap_zip', default='gz',
    type=click.Choice(['none', 'gz', 'bz2', 'xz'], case_sensitive=False),
    help='Compress read-to-feature maps using this algorithm.')
# performance
@click.option(
    '--chunk', type=click.INT, default=None,
    help='Number of alignment lines to read and parse in each chunk.')
@click.option(
    '--cache', type=click.INT, default=1024,
    help='Number of recent results to cache for faster classification.')
@click.option(
    '--no-exe', is_flag=True,
    help='Disable calling external programs for decompression.')
def classify_cmd(**kwargs):
    """Generate a profile of samples based on a classification system.
    """
    workflow(**kwargs)


# `tools` provides utilities for working with alignments, maps and profiles

@cli.group('tools', cls=NaturalOrderGroup)
def tools():
    """Utilities for working with alignments, maps and profiles.
    """
    pass  # pragma: no cover


@tools.command('filter')
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to input profile.')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output profile.')
@click.option(
    '--min-count', '-c', type=click.IntRange(min=1),
    help='Per-sample minimum count threshold.')
@click.option(
    '--min-percent', '-p', type=click.FLOAT,
    help='Per-sample minimum percentage threshold.')
@click.pass_context
def filter_cmd(ctx, **kwargs):
    """Filter a profile by per-sample abundance.
    """
    filter_wf(**kwargs)


@tools.command('merge')
@click.option(
    '--input', '-i', 'input_fps', required=True, multiple=True,
    type=click.Path(exists=True),
    help=('Path to input profiles or directories containing profiles. Can '
          'accept multiple paths.'))
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output profile.')
@click.pass_context
def merge_cmd(ctx, **kwargs):
    """Merge multiple profiles into one profile.
    """
    merge_wf(**kwargs)


@tools.command('collapse')
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to input profile.')
@click.option(
    '--map', '-m', 'map_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help=('Mapping of source features to target features. (supports '
          'many-to-many relationships).'))
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output profile.')
@click.option(
    '--normalize', '-z', is_flag=True,
    help=('Count each target feature as 1/k (k is the number of targets '
          'mapped to a source). Otherwise, count as one.'))
@click.option(
    '--names', '-n', 'names_fp', type=click.Path(exists=True),
    help='Names of target features to append to the output profile.')
@click.pass_context
def collapse_cmd(ctx, **kwargs):
    """Collapse a profile based on feature mapping.
    """
    collapse_wf(**kwargs)


@tools.command('coverage')
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to input profile.')
@click.option(
    '--map', '-m', 'map_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Mapping of feature groups to member features.')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output coverage table.')
@click.option(
    '--threshold', '-t', type=click.IntRange(1, 100),
    help=('Convert coverage to presence (1) / absence (0) data by this '
          'percentage threshold.'))
@click.option(
    '--count', '-c', is_flag=True,
    help=('Record numbers of covered features instead of percentages '
          '(overrides threshold).'))
@click.option(
    '--names', '-n', 'names_fp', type=click.Path(exists=True),
    help='Names of feature groups to append to the coverage table.')
@click.pass_context
def coverage_cmd(ctx, **kwargs):
    """Calculate per-sample coverage of feature groups.
    """
    coverage_wf(**kwargs)


if __name__ == '__main__':
    cli()  # pragma: no cover
