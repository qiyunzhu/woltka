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


def profile_to_biom(profile, samples=None, tree=None, rankd=None, named=None,
                    add_lineage=False):
    """Convert a profile into a BIOM table.

    Parameters
    ----------
    profile : dict
        Profile to convert.
    samples : list of str, optional
        Sample ID list to include.
    tree : dict, optional
        Taxonomic tree.
    rankd : dict, optional
        Rank dictionary.
    named : dict, optional
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
            'Name': named[feature] if named and feature in named else None,
            'Rank': rankd[feature] if rankd and feature in rankd else None,
            'Lineage': get_lineage(feature, tree) if add_lineage else None})
    return biom.Table(np.array(data), features, samples, metadata,
                      generated_by=f'woltka-{__version__}')


def write_biom(table: biom.Table, path_: str):
    """Write a BIOM table to file.

    Parameters
    ----------
    table : biom.Table
        BIOM table to write.
    path_ : str
        Output filepath.
    """
    with biom.util.biom_open(path_, 'w') as f:
        table.to_hdf5(f, table.generated_by)
