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
from .tree import get_lineage_gg
from .__init__ import __version__


def profile_to_biom(profile, samples=None, tree=None, rankdic=None,
                    namedic=None, name_as_id=False):
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
    name_as_id : bool, optional
        Replace feature IDs with names. It applies to row headers and "Lineage"
        column, and removes "Name" column.

    Returns
    -------
    biom.Table
        Converted BIOM table.
    """
    features = sorted(allkeys(profile))
    samples = samples or sorted(profile)
    data, metadata, names = [], [], []
    for feature in features:

        # generate data
        row = []
        for sample in samples:
            row.append(str(profile[sample][feature]) if feature in
                       profile[sample] else 0)
        data.append(row)

        # generate metadata
        meta_ = {}
        if namedic:
            name = namedic[feature] if feature in namedic else None
            if name_as_id:
                names.append(name or feature)
            else:
                meta_['Name'] = name or None
        if rankdic:
            meta_['Rank'] = rankdic[feature] if feature in rankdic else None
        if tree:
            meta_['Lineage'] = get_lineage_gg(
                feature, tree, namedic if name_as_id else None) or None
        if meta_:
            metadata.append(meta_)

    # build biom table
    return biom.Table(np.array(data), names or features, samples, metadata or
                      None, generated_by=f'woltka-{__version__}')


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
