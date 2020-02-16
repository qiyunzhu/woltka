#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import join
import click

from woltka.core import count, assign
from woltka.util import readzip, id2file_map, allkeys
from woltka.parse import read_map_file
from woltka.tree import build_tree


@click.command()
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help='directory of input read alignment(s)')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True),
    help=('path to output profile file (single rank) or directory (multiple '
          'ranks)'))
# input information
@click.option(
    '--format', '-f', 'input_fmt', default='auto',
    type=click.Choice(['auto', 'b6o', 'sam', 'map'], case_sensitive=False),
    help=('format of read alignment: "auto": automatic determination '
          '(default), "b6o": BLAST tabular format (-outfmt 6), "sam": SAM '
          'format, "map": simple map of query <tab> subject'))
@click.option(
    '--extension', '-e', 'input_ext',
    help='input filename extension following sample ID')
@click.option(
    '--sample-ids', '-s', type=click.File('r'),
    help='list of sample IDs to be included')
# behavior
@click.option(
    '--rank', '-r', 'rank_lst', type=click.STRING,
    help=('classify sequences at this rank; ignore or enter "none" to omit '
          'classification; enter "free" for free-rank classification; can '
          'specify multiple comma-delimited ranks and one profile will be '
          'generated for each rank'))
@click.option(
    '--multi/--no-multi', default=True,
    help=('allow one sequence to be assigned to multiple classification '
          'units at the same rank; per-unit match counts will be recorded '
          'and profile will be normalized by total number of matches'))
@click.option(
    '--lca/--no-lca', default=True,
    help=('attempt to find lowest common ancestor (LCA) in taxonomy for '
          'non-unique matches; note: root is not assumed unless defined, '
          'and there may be multiple LCAs, in which case each is counted '
          'once; results are then subject to ambiguity treatment'))
@click.option(
    '--ixend/--no-ixend', default=False,
    help=('subject identifiers end with underscore index, the latter of which '
          'is to be removed prior to mapping.'))
# tree information
@click.option(
    '--map', '-m', 'map_fps', type=click.Path(exists=True), multiple=True,
    help=('map(s) of subjects to higher classification units, such as '
          'nucleotides to host genomes, sequence IDs to taxonomy IDs, gene '
          'family to pathway, etc., can accept multiple maps entered in low-'
          'to-high order'))
@click.option(
    '--names', 'names_fp', type=click.Path(exists=True),
    help=('map of taxonomic units to labels; can be plain map or NCBI '
          'names.dmp'))
@click.option(
    '--nodes', 'nodes_fp', type=click.Path(exists=True),
    help=('hierarchical structure of taxonomy defined by NCBI names.dmp '
          'or compatible format'))
@click.option(
    '--newick', 'newick_fp', type=click.Path(exists=True),
    help=('classification hierarchies defined by a tree in Newick format.'))
@click.option(
    '--rank-table', 'ranktb_fp', type=click.Path(exists=True),
    help=('classification hierarchies defined by a table with each column '
          'representing a level and header as level name.'))
@click.option(
    '--lineage', 'lineage_fp', type=click.Path(exists=True),
    help=('map of subjects/groups to lineage strings in format of "taxonomic;'
          'units;from;high;to;low", can be Greengenes-style taxonomy where '
          'level codes such as "k__" will be parsed'))
def classify(input_fp, output_fp,
             input_fmt, input_ext, sample_ids,
             rank_lst, multi, lca, ixend,
             map_fps, names_fp, nodes_fp, newick_fp, ranktb_fp, lineage_fp):
    """Generate a profile of query samples based on hierarchical organization
    subjects."""

    # parse sample Ids
    samples = None
    if sample_ids:
        samples = sample_ids.read().splitlines()
        click.echo('Samples to include: %d.' % len(samples))

    # match input files with sample Ids
    input_map = id2file_map(input_fp, input_ext, samples)
    if not samples:
        samples = sorted(input_map.keys())
        click.echo('Samples to read: %d.' % len(samples))
    elif len(input_map) < len(samples):
        raise ValueError('Inconsistent sample IDs.')

    # build classification system
    tree, rankd, named, root = build_tree(
        map_fps, names_fp, nodes_fp, newick_fp, ranktb_fp, lineage_fp)

    # parse target ranks
    ranks = ['none'] if rank_lst is None else rank_lst.split(',')
    data = {x: {} for x in ranks}

    # parse input maps and generate profile
    args = [tree, rankd, root]
    for sample in samples:
        click.echo(f'Parsing {sample}...')

        # read alignment file into query-subject(s) map
        with readzip(join(input_fp, input_map[sample])) as f:
            map_ = read_map_file(f, input_fmt)

        # merge duplicate query-subject pairs
        map_ = {k: set(v) for k, v in map_.items()}
        click.echo(f'Query sequences: {len(map_)}.')
        for rank in ranks:
            assignment = {k: assign(v, rank, *args) for k, v in map_.items()}
            data[rank][sample] = count(assignment)
            try:
                del data[rank][sample][None]
            except KeyError:
                pass
        # click.echo(data[sample].keys())

    # determine output filenames
    if len(ranks) == 1:
        rank2fp = {ranks[0]: output_fp}
    else:
        makedirs(output_fp, exist_ok=True)
        rank2fp = {x: join(output_fp, f'{x}.tsv') for x in ranks}

    # write output profile(s)
    for rank in ranks:
        with open(rank2fp[rank], 'w') as f:
            write_profile(f, data[rank], named, samples)
    click.echo('Done.')


def write_profile(fh, data, named=None, samples=None):
    """Write profile to a plain tab-delimited file.
    """
    if samples is None:
        samples = sorted(data)
    print('#SampleID\t{}'.format('\t'.join(samples)), file=fh)
    for key in sorted(allkeys(data)):
        # get feature name
        try:
            row = [named[key]]
        except (TypeError, KeyError):
            row = [key]
        # get feature count
        for sample in samples:
            try:
                row.append(str(int(data[sample][key])))
            except KeyError:
                row.append('0')
        # skip all-zero row
        if any(x != '0' for x in row[1:]):
            print('\t'.join(row), file=fh)


if __name__ == "__main__":
    classify()
