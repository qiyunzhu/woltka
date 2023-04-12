#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Features for dealing with ambiguous matches.
"""

from sys import exit
from os.path import join
import numpy as np
import click
from biom import Table

from .table import read_table, table_shape, write_table
from .file import id2file_from_dir, readzip


def disambig_wf(input_fp:    str,
                map_dir:     str,
                output_fp:   str,
                min_count:   int = None,
                min_ratio: float = None):
    """Workflow for dropping features that have inadequate unique matches from
    sequencing reads.

    Raises
    ------
    SystemExit
        Both thresholds are provided.

    See Also
    --------
    .cli.filter_cmd
        Command-line arguments and help information.
    """
    # validate parameters
    if all((min_count, min_ratio)):
        exit('Only one of minimum count or minimum ratio thresholds '
             'can be specified.')

    # read input profile
    table, _ = read_table(input_fp)
    n = table_shape(table)[0]
    click.echo(f'Number of features before disambiguation: {n}.')

    # get sample list
    samples = table.ids().tolist() if isinstance(table, Table) else table[2]

    # validate read mappings
    mapfiles = id2file_from_dir(map_dir, ids=samples)

    # available external compressors
    zippers = {}

    # filter profile by threshold
    click.echo('Disambiguate profile...')
    todels = []
    for sample in samples:
        fname = mapfiles[sample]
        features = {}

        # read mappings
        with readzip(join(map_dir, fname), zippers) as f:
            for line in f:
                _, _, targets = line.rstrip().partition('\t')

                # ambiguous matches
                if '\t' in targets:
                    for target in targets.split('\t'):
                        target = target.rsplit(':', 1)[0]
                        features.setdefault(target, [0, 0])[1] += 1

                # unique match
                else:
                    features.setdefault(targets, [0, 0])[0] += 1

        # drop features
        todel = []
        todel_append = todel.append
        if min_count:
            for feature, [uniq, _] in features.items():
                if uniq < min_count:
                    todel_append(feature)
        elif min_ratio:
            for feature, [uniq, ambi] in features.items():
                if uniq / (ambi + uniq) < min_ratio:
                    todel_append(feature)
        else:
            for feature, [_, ambi] in features.items():
                if ambi:
                    todel_append(feature)
        todels.append(set(todel))
        click.echo(f'{sample}: Dropped {len(todel)} out of {len(features)}'
                   ' features.')

    # filter table
    if isinstance(table, Table):
        features = table.ids('observation').tolist()
        metadata = table.metadata(axis='observation')
        data = table.matrix_data.toarray()
        for i, feature in enumerate(features):
            mask = np.array([feature in x for x in todels])
            data[i, mask] = 0
        table = Table(data, features, samples, metadata)
        table.remove_empty(axis='observation')
    else:
        empty = []
        for i, feature in enumerate(table[1]):
            row = table[0][i]
            for j, todel in enumerate(todels):
                if row[j] and feature in todel:
                    row[j] = 0
            if not any(row):
                empty.append(i)
        for i in reversed(empty):
            del table[0][i]
            del table[1][i]
            del table[3][i]

    click.echo(' Done.')
    n = table_shape(table)[0]
    click.echo(f'Number of features after disambiguation: {n}.')

    # write disambiguated profile
    write_table(table, output_fp)
    click.echo('Disambiguated profile written.')
