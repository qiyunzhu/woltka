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

from woltk.core import count
from woltk.util import readzip, id2file_map, allkeys
from woltk.parse import read_map_file


@click.command()
@click.option(
    '--input-dir', '-i', required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help='input read map directory')
@click.option(
    '--output-table', '-o', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='output profile')
@click.option(
    '--input-fmt', '-f', default='auto',
    type=click.Choice(['auto', 'b6o', 'sam', 'map'], case_sensitive=False),
    help=('format of read alignment: "auto": automatic determination '
          '(default), "b6o": BLAST tabular format (-outfmt 6), "sam": SAM '
          'format, "map": simple map of query <tab> subject'))
@click.option(
    '--input-ext', '-e',
    help='input filename extension following sample ID')
@click.option(
    '--sample-ids', '-s', type=click.File('r'),
    help='list of sample IDs to be included')
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
@click.option(
    '--groups', '-g', type=click.Path(exists=True), multiple=True,
    help=('map(s) of subjects to higher groups, such as nucleotides to host '
          'genomes, or sequence IDs to taxonomy IDs, can accept multiple maps '
          'specified in order'))
@click.option(
    '--lineage', '-l',
    type=click.Path(writable=True, dir_okay=False),
    help=('map of subjects/groups to lineage strings in format of "taxonomic;'
          'units;from;high;to;low", can be Greengenes-style taxonomy where '
          'level codes such as "k__" will be parsed'))
@click.option(
    '--nodes', type=click.Path(exists=True),
    help=('hierarchical structure of taxonomy defined by NCBI names.dmp '
          'or compatible format'))
@click.option(
    '--names', type=click.Path(exists=True),
    help=('map of taxonomic units to labels; can be plain map or NCBI '
          'names.dmp'))
def classify(input_dir, output_table, input_fmt, input_ext, sample_ids,
             ambiguity, lca, ixend, groups, lineage, nodes, names):
    """Generate a profile of query samples based on hierarchical organization
    subjects."""

    # parse sample Ids
    ids = None
    if sample_ids:
        ids = sample_ids.read().splitlines()
        click.echo('Samples to include: %d.' % len(ids))

    # match input files with sample Ids
    input_map = id2file_map(input_dir, input_ext, ids)
    if not ids:
        ids = sorted(input_map.keys())
        click.echo('Samples to read: %d.' % len(ids))
    elif len(input_map) < len(ids):
        raise ValueError('Inconsistent sample IDs.')

    # read maps of higher groups
    group_maps = []
    for group in groups:
        group_maps.append(readzip(group))

    # parse input maps and generate profile
    data = {}
    for id_ in ids:
        click.echo(f'Parsing {id_}...')
        with readzip(join(input_dir, input_map[id_])) as f:
            x = read_map_file(f, input_fmt)
            click.echo(len(x))
            data[id_] = count(x, ambiguity)
            click.echo(data[id_].keys())

    # write output profile
    with open(output_table, 'w') as f:
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
