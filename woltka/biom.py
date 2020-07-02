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
from .tree import lineage_str
from .__init__ import __name__, __version__


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
    if samples:
        samples = [x for x in samples if x in profile]
    else:
        samples = sorted(profile)
    observations, data, metadata = [], [], []

    for key in sorted(allkeys(profile)):
        # generate data
        row = []
        for sample in samples:
            row.append(profile[sample][key] if key in profile[sample] else 0)
        data.append(row)
        # generate metadata and observation
        stratum, feature = key if isinstance(key, tuple) else (None, key)
        obsrv, meta_ = feature, {}
        # get feature name
        if namedic:
            name = namedic[feature] if feature in namedic else None
            if name_as_id:
                obsrv = name or feature
            else:
                meta_['Name'] = name or None
        # prefix with stratum label
        if stratum:
            obsrv = f'{stratum}|{obsrv}'
        observations.append(obsrv)
        # get feature rank
        if rankdic:
            meta_['Rank'] = rankdic[feature] if feature in rankdic else None
        # get lineage string
        if tree:
            meta_['Lineage'] = lineage_str(
                feature, tree, namedic if name_as_id else None) or None
        if meta_:
            metadata.append(meta_)

    # build biom table
    return biom.Table(np.array(data), observations, samples, metadata or None,
                      generated_by=f'{__name__}-{__version__}')


def table_to_biom(data, observations, samples, metadata=None):
    """Convert a profile into a BIOM table.

    Parameters
    ----------
    data : list of list
        Table data.
    observations : list
        Observation IDs.
    samples : list
        Sample IDs.
    metadata : list of dict, optional
        Observation metadata.

    Returns
    -------
    biom.Table
        Converted BIOM table.
    """
    return biom.Table(np.array(data), observations, samples, metadata or None,
                      generated_by=f'{__name__}-{__version__}')


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


def filter_biom(table: biom.Table, th: float):
    """Filter a BIOM table by per-sample percentage frequency threshold.

    Parameters
    ----------
    table : biom.Table
        BIOM table to filter.
    th : float
        Per-sample minimum abundance threshold. If >= 1, this value is an
        absolute count; if < 1, it is a fraction of sum of counts.

    Returns
    -------
    biom.Table
        Filtered BIOM table.
    """
    def f(data, id_, md):
        bound = th if th > 1 else data.sum() * th
        data[data < bound] = 0
        return data

    res = table.copy()
    res.transform(f, axis='sample')
    res.remove_empty(axis='observation')
    return res
