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
from .tools import normalize_wf, filter_wf, merge_wf, collapse_wf, coverage_wf


class NaturalOrderGroup(click.Group):
    """Natural ordering of click command groups.
    """
    def list_commands(self, ctx):
        return self.commands.keys()


GRP_KA = dict(
    cls=NaturalOrderGroup,
    context_settings=dict(help_option_names=['-h', '--help']))
CMD_KA = dict(
    no_args_is_help=True)


@click.version_option(__version__)
@click.group(**GRP_KA)
def cli():
    """Woltka: a versatile meta'omic data classifier.
    """
    pass  # pragma: no cover


@cli.command('classify', **CMD_KA)
# input and output
@click.option(
    '--input', '-i', 'input_fp', required=True, type=click.Path(
        exists=True, file_okay=True, dir_okay=True, allow_dash=True),
    help=('Path to input alignment file or directory of alignment files.'
          ' Enter "-" for stdin.'))
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
    '--trim-sub', 'trimsub', is_flag=False, flag_value='_',
    help=('Trim subject IDs at the last given delimiter. Default: "_", or '
          'enter a custom value.'))
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
    '--map-as-rank/--map-no-rank', 'map_rank', default=None,
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
# normalization
@click.option(
    '--sizes', '-z', type=click.Path(exists=True),
    help=('Divide counts by subject sizes. Can provide a mapping file, or '
          'type "." to calculate from gene coordinates.'))
@click.option(
    '--frac', is_flag=True,
    help='Divide counts by total count of each sample (i.e., fractions).')
@click.option(
    '--scale', type=click.STRING,
    help='Scale counts by this factor. Accepts "k", "M" suffixes.')
@click.option(
    '--digits', type=click.IntRange(0, 10),
    help='Round counts to this number of digits after the decimal point.')
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
@click.option(
    '--outcov', 'outcov_dir', type=click.Path(dir_okay=True),
    help='Write subject coverage maps to this directory.')
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
    """Main classification workflow: Alignments => profile(s).
    """
    workflow(**kwargs)


@cli.command('collapse', **CMD_KA)
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to input profile.')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output profile.')
@click.option(
    '--map', '-m', 'map_fp',
    type=click.Path(exists=True, dir_okay=False),
    help=('Mapping of source features to target features. Supports '
          'many-to-many relationships.'))
@click.option(
    '--divide', '-d', is_flag=True,
    help=('Count each target feature as 1/k (k is the number of targets '
          'mapped to a source). Otherwise, count as one.'))
@click.option(
    '--field', '-f', type=click.INT,
    help=('Collapse x-th field of stratified features. For example, "A|a" '
          'has fields 1 ("A") and 2 ("a").'))
@click.option(
    '--nested', '-e', is_flag=True,
    help=('Fields are nested (each field is a child of the previous field). '
          'For example, "A_1" represents "1" of "A".'))
@click.option(
    '--sep', '-s', type=click.STRING,
    help=('Field separator for nested features (default: "_") or otherwise '
          '(default: "|").'))
@click.option(
    '--names', '-n', 'names_fp', type=click.Path(exists=True),
    help='Names of target features to append to the output profile.')
def collapse_cmd(**kwargs):
    """Collapse a profile by feature mapping and/or hierarchy.
    """
    collapse_wf(**kwargs)


@cli.command('normalize', **CMD_KA)
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to input profile.')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output profile.')
@click.option(
    '--sizes', '-z', 'sizes_fp',
    type=click.Path(exists=True, dir_okay=False),
    help=('Path to mapping of feature sizes, by which values will be divided. '
          'If omitted, will divide values by sum per sample.'))
@click.option(
    '--scale', '-s', type=click.STRING,
    help='Scale values by this factor. Accepts "k", "M" suffixes.')
@click.option(
    '--digits', '-d', type=click.IntRange(0, 10),
    help=('Round values to this number of digits after the decimal point. If '
          'omitted, will keep decimal precision of input profile.'))
def normalize_cmd(**kwargs):
    """Normalize a profile to fractions and/or by feature sizes.
    """
    normalize_wf(**kwargs)


@cli.command('filter', **CMD_KA)
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
def filter_cmd(**kwargs):
    """Filter a profile by per-sample abundance.
    """
    filter_wf(**kwargs)


@cli.command('merge', **CMD_KA)
@click.option(
    '--input', '-i', 'input_fps', required=True, multiple=True,
    type=click.Path(exists=True),
    help=('Path to input profiles or directories containing profiles. Can '
          'accept multiple paths.'))
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='Path to output profile.')
def merge_cmd(**kwargs):
    """Merge multiple profiles into one profile.
    """
    merge_wf(**kwargs)


@cli.command('coverage', **CMD_KA)
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
def coverage_cmd(**kwargs):
    """Calculate per-sample coverage of feature groups.
    """
    coverage_wf(**kwargs)


# the "tools" menu is for backward compatibility
@cli.group('tools', cls=NaturalOrderGroup)
def tools():
    """Entries to the same commands for backward compatibility (deprecated).
    """
    pass  # pragma: no cover


for cmd in collapse_cmd, normalize_cmd, filter_cmd, merge_cmd, coverage_cmd:
    tools.add_command(cmd)


if __name__ == '__main__':
    cli()  # pragma: no cover
