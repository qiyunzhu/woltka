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
from .__init__ import __name__, __version__


def table_to_biom(data, observations, samples, metadata=None):
    """Convert table components into a BIOM table.

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
    return biom.Table(np.array(data), observations, samples, metadata or None)


def biom_to_table(table: biom.Table):
    """Convert a BIOM table into table components.

    Parameters
    ----------
    table : biom.Table
        BIOM table to convert.

    Returns
    -------
    tuple of
        list of list
            Profile data.
        list
            Observation IDs.
        list
            Sample IDs.
        list of dict
            Observation metadata.
    """
    return (table.matrix_data.toarray().astype(int).tolist(),
            table.ids('observation').tolist(),
            table.ids('sample').tolist(),
            list(map(dict, table.metadata(axis='observation') or ())))


def write_biom(table: biom.Table, fp: str):
    """Write a BIOM table to file.

    Parameters
    ----------
    table : biom.Table
        BIOM table to write.
    fp : str
        Output filepath.

    Notes
    -----
    The `generated_by` attribute of the output BIOM table will be like
    "woltka-version".
    """
    with biom.util.biom_open(fp, 'w') as f:
        table.to_hdf5(f, f'{__name__}-{__version__}')


def filter_biom(table: biom.Table, th: float):
    """Filter a BIOM table by per-sample count or percentage threshold.

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
        bound = th if th >= 1 else data.sum() * th
        data[data < bound] = 0
        return data

    res = table.copy()
    res.transform(f, axis='sample')
    res.remove_empty(axis='observation')
    return res
