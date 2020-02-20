#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from . import __version__
from .classify import classify as _classify


class NatOrder(click.Group):
    """Natural ordering of click command groups.

    See Also
    --------
    https://github.com/pallets/click/issues/513
    """
    def list_commands(self, ctx):
        return self.commands.keys()


@click.version_option(__version__)
@click.group(cls=NatOrder)
def cli():
    pass


# `gotu` is a simplified wrapper of `classify`

@cli.command('gotu')
@click.option(
    '--input', '-i', 'input_path', required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    help=('Path to a multiplexed alignment file or a directory of per-sample '
          'alignment files.'))
@click.option(
    '--output', '-o', 'output_path', required=True,
    type=click.Path(writable=True),
    help=('Path to output gOTU table.'))
@click.option(
    '--ambig/--no-ambig', default=True,
    help=('Allow one sequence to be assigned to multiple gOTUs. Each hit '
          'will be counted as 1 / k (k is the totally number of hits).'))
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
    '--input', '-i', 'input_path', required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    help='Path to input alignment file or directory of alignment files.')
@click.option(
    '--output', '-o', 'output_path', required=True,
    type=click.Path(writable=True),
    help='Path to output profile file or directory of profile files.')
# input information
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
    '--sample-ids', '-s', type=click.File('r'),
    help='List of sample IDs to be included.')
@click.option(
    '--demux/--no-demux', default=None,
    help='Demultiplex alignment by first underscore in query identifier.')
# classification
@click.option(
    '--rank', '-r', 'rank_lst', type=click.STRING,
    help=('Classify sequences at this rank. Ignore or enter "none" to omit '
          'classification; enter "free" for free-rank classification. Can '
          'specify multiple comma-delimited ranks and one profile will be '
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
    help='Can report subject IDs in classification result.')
@click.option(
    '--lca/--no-lca', default=True,
    help='Find lowest common ancestor (LCA) for non-unique matches.')
@click.option(
    '--deidx/--no-deidx', default=False,
    help='Strip "underscore index" suffixes from subject IDs.')
# gene information
@click.option(
    '--coords', '-c', 'coords_fp', type=click.Path(exists=True),
    help='Table of gene coordinates of  on reference genomes.')
# tree information
@click.option(
    '--map', '-m', 'map_fps', type=click.Path(exists=True), multiple=True,
    help=('Map(s) of subjects or lower classification units to higher ones. '
          'Can accept multiple maps.'))
@click.option(
    '--nodes', 'nodes_fp', type=click.Path(exists=True),
    help='Hierarchies defined by NCBI nodes.dmp or compatible formats.')
@click.option(
    '--newick', 'newick_fp', type=click.Path(exists=True),
    help='Hierarchies defined by a tree in Newick format.')
@click.option(
    '--rank-table', 'ranktb_fp', type=click.Path(exists=True),
    help='Table of classification units at each rank (column).')
@click.option(
    '--lineage', 'lineage_fp', type=click.Path(exists=True),
    help='Map of lineage strings. Can accept Greengenes-style rank prefix.')
@click.option(
    '--names', 'names_fp', type=click.Path(exists=True),
    help=('Names of classification units as defined by NCBI names.dmp or a '
          'plain map.'))
def classify(**kwargs):
    """Generate a profile of samples based on a classification system.
    """
    _classify(**kwargs)


if __name__ == '__main__':
    cli()
