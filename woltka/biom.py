#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for operating BIOM tables.
"""

import numpy as np
import biom
from .util import allkeys
from .tree import get_lineage
from .__init__ import __version__


def profile_to_biom(profile, samples=None, tree=None, rankdic=None,
                    namedic=None, add_lineage=False):
    """Convert a profile into a BIOM table.

    Parameters
    ----------
    profile : dict
        Profile to convert.
    samples : list of str, optional
        Sample ID list to include.
    tree : dict, optional
        Taxonomic tree.
    rankdic : dict, optional
        Rank dictionary.
    namedic : dict, optional
        Taxon name dictionary.
    add_lineage : bool, optional
        Append lineage (root-to-current hierarchies) to metadata.

    Returns
    -------
    biom.Table
        Converted BIOM table.
    """
    features = sorted(allkeys(profile))
    samples = samples or sorted(profile)
    data, metadata = [], []
    for feature in features:
        row = []
        for sample in samples:
            try:
                row.append(profile[sample][feature])
            except KeyError:
                row.append(0)
        data.append(row)
        metadata.append({
            'Name': (
                namedic[feature] if namedic and feature in namedic else None),
            'Rank': (
                rankdic[feature] if rankdic and feature in rankdic else None),
            'Lineage': (
                get_lineage(feature, tree) if add_lineage else None)})
    return biom.Table(np.array(data), features, samples, metadata,
                      generated_by=f'woltka-{__version__}')


def write_biom(table: biom.Table, fp: str):
    """Write a BIOM table to file.

    Parameters
    ----------
    table : biom.Table
        BIOM table to write.
    fp : str
        Output filepath.
    """
    with biom.util.biom_open(fp, 'w') as f:
        table.to_hdf5(f, table.generated_by)
