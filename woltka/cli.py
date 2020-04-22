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


class NaturalOrderGroup(click.Group):
    """Natural ordering of click command groups.
    """
    def list_commands(self, ctx):
        return self.commands.keys()


@click.version_option(__version__)
@click.group(cls=NaturalOrderGroup)
def cli():
    pass


# `gotu` is a simplified wrapper of `classify`

@cli.command('gotu')
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    help=('Path to a multiplexed alignment file or a directory of per-sample '
          'alignment files.'))
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True),
    help=('Path to output gOTU table.'))
@click.option(
    '--ambig/--no-ambig', default=True,
    help=('Allow one sequence to be assigned to multiple gOTUs. Each hit '
          'will be counted as 1 / k (k is the totally number of hits). '
          'Otherwise, sequences with multiple matches will be dropped.'))
@click.option(
    '--map', '-m', 'map_fps', type=click.Path(exists=True), multiple=True,
    help=('Map of nucleotides to genomes.'))
@click.pass_context
def gotu(ctx, **kwargs):
    """Generate a gOTU table based on sequence alignments.
    """
    ctx.invoke(classify, **kwargs)


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
    '--lines', type=click.INT, default=1000000, show_default=True,
    help=('Number of lines to read from alignment file per chunk.'))
# hierarchies
@click.option(
    '--nodes', 'nodes_fp', type=click.Path(exists=True),
    help='Hierarchies defined by NCBI nodes.dmp or compatible formats.')
@click.option(
    '--newick', 'newick_fp', type=click.Path(exists=True),
    help='Hierarchies defined by a tree in Newick format.')
@click.option(
    '--lineage', 'lineage_fp', type=click.Path(exists=True),
    help='Map of lineage strings. Can accept Greengenes-style rank prefix.')
@click.option(
    '--rank-table', 'rank_table_fp', type=click.Path(exists=True),
    help='Table of classification units at each rank (column).')
@click.option(
    '--map', '-m', 'map_fps', type=click.Path(exists=True), multiple=True,
    help=('Map of lower classification units to higher ones. Can accept '
          'multiple files.'))
@click.option(
    '--map-as-rank', is_flag=True,
    help='Map filename stem is rank name.')
@click.option(
    '--names', 'names_fps', type=click.Path(exists=True), multiple=True,
    help=('Names of classification units as defined by NCBI names.dmp or a '
          'simple map. Can accept multiple files.'))
# assignment
@click.option(
    '--rank', '-r', 'ranks', type=click.STRING,
    help=('Classify sequences at this rank. Ignore or enter "none" to omit '
          'classification; enter "free" for free-rank classification. Can '
          'specify multiple comma-separated ranks and one profile will be '
          'generated for each rank.'))
@click.option(
    '--above/--no-above', default=False,
    help='Allow assigning to a classification unit higher than given rank.')
@click.option(
    '--major', type=click.IntRange(51, 99),
    help='Majority-rule assignment percentage threshold.')
@click.option(
    '--ambig/--no-ambig', default=True,
    help='Allow assigning one sequence to multiple classification units.')
@click.option(
    '--subok/--no-subok', default=True,
    help='Allow assigning sequences to their subjects.')
@click.option(
    '--deidx/--no-deidx', default=False,
    help='Strip "_index" suffixes from subject IDs.')
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
    help='Output feature table format (BIOM or TSV).')
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
    '--outmap-zip', default='gz',
    type=click.Choice(['none', 'gz', 'bz2', 'xz'], case_sensitive=False),
    help=('Compress read maps using this algorithm.'))
def classify(**kwargs):
    """Generate a profile of samples based on a classification system.

    Notes
    -----
    Details of parameters are provided in `workflow.py` and `doc/cli.md`.

    See Also
    --------
    workflow.workflow
    """
    workflow(**kwargs)


if __name__ == '__main__':
    cli()
