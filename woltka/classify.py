#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
import click

from woltka.core import count
from woltka.util import readzip, id2file_map, allkeys
from woltka.parse import read_map_file
from woltka.tree import build_tree


@click.command()
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help='input read map directory')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='output profile')
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
    '--ambiguity', '-a', default='norm', show_default=True,
    type=click.Choice(['uniq', 'all', 'norm'], case_sensitive=False),
    help=('how to treat reads that occur multiple times: "uniq": drop '
          'non-unique reads; "all": count each occurrence once; "norm": '
          'count each occurrence 1/k times (k = number of occurrences)'))
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
    '--levels', 'levels_fp', type=click.Path(exists=True),
    help=('classification hierarchies defined by a table with each column '
          'representing a level and header as level name.'))
@click.option(
    '--lineage', 'lineage_fp', type=click.Path(exists=True),
    help=('map of subjects/groups to lineage strings in format of "taxonomic;'
          'units;from;high;to;low", can be Greengenes-style taxonomy where '
          'level codes such as "k__" will be parsed'))
@click.option(
    '--groups', 'group_fps', type=click.Path(exists=True), multiple=True,
    help=('map(s) of subjects to higher groups, such as nucleotides to host '
          'genomes, or sequence IDs to taxonomy IDs, can accept multiple maps '
          'specified in order'))
def classify(input_fp, output_fp,
             input_fmt, input_ext, sample_ids,
             ambiguity, lca, ixend,
             names_fp, nodes_fp, newick_fp, levels_fp, lineage_fp, group_fps):
    """Generate a profile of query samples based on hierarchical organization
    subjects."""

    '''Read sample information.'''

    # parse sample Ids
    ids = None
    if sample_ids:
        ids = sample_ids.read().splitlines()
        click.echo('Samples to include: %d.' % len(ids))

    # match input files with sample Ids
    input_map = id2file_map(input_fp, input_ext, ids)
    if not ids:
        ids = sorted(input_map.keys())
        click.echo('Samples to read: %d.' % len(ids))
    elif len(input_map) < len(ids):
        raise ValueError('Inconsistent sample IDs.')

    '''Read classification tree.'''

    tree, ranks, names = build_tree(
        names_fp, nodes_fp, newick_fp, levels_fp, lineage_fp, group_fps)

    '''Parse maps.'''

    # parse input maps and generate profile
    data = {}
    for id_ in ids:
        click.echo(f'Parsing {id_}...')
        with readzip(join(input_fp, input_map[id_])) as f:
            x = read_map_file(f, input_fmt)
            click.echo(len(x))
            data[id_] = count(x, ambiguity)
            click.echo(data[id_].keys())

    # write output profile
    with open(output_fp, 'w') as f:
        write_profile(f, data, samples=ids)
    click.echo('Done.')


def write_profile(fh, data, fmt='tsv', samples=None):
    """Write profile to a plain tab-delimited file."""
    if samples is None:
        samples = sorted(data)
    if fmt == 'tsv':
        write_tsv(fh, data, samples)
    elif fmt == 'biom':
        write_biom(fh, data, samples)


def write_tsv(fh, data, samples):
    print('#SampleID\t{}'.format('\t'.join(samples)), file=fh)
    for feature in sorted(allkeys(data)):
        row = [feature]
        for sample in samples:
            try:
                row.append(str(int(data[sample][feature])))
            except KeyError:
                row.append('0')
        print('\t'.join(row), file=fh)


def write_biom(fh, data, samples):
    pass


if __name__ == "__main__":
    classify()
