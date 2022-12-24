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

from itertools import accumulate
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


def biom_max_f(table: biom.Table):
    """Get the maximum number of digits after the decimal place in a table.

    Parameters
    ----------
    table : biom.Table
        BIOM table to examine.

    Returns
    -------
    int
        Maximum decimal precision.
    """
    return max([0] + [str(x)[::-1].find('.') for x in table.matrix_data.data])


def divide_biom(table: biom.Table, sizes: dict):
    """Divide table values by feature sizes in place.

    Parameters
    ----------
    table : biom.Table
        BIOM table to divide.
    sizes : dict
        Feature size dictionary.
    """
    table.transform(lambda data, id_, md: data / sizes[id_],
                    axis='observation', inplace=True)


def scale_biom(table: biom.Table, scale: float):
    """Multiply table values by a constant factor in place.

    Parameters
    ----------
    table : biom.Table
        BIOM table to multiply.
    scale : float
        Scale factor.
    """
    table.transform(lambda data, id_, md: data * scale, inplace=True)


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


def round_biom(table: biom.Table, digits=0):
    """Round a BIOM table's data and drop empty observations in place.

    Parameters
    ----------
    table : biom.Table
        BIOM table to round.
    digits : int, optional
        Digits after the decimal point.

    Notes
    -----
    There is a fully vectorized, much faster alternate:

    >>> arr = table.matrix_data.data
    >>> near = (arr * 2).round(digits) / 2
    >>> choice = np.abs(arr - near) <= error
    >>> table.matrix_data.data = np.where(choice, near, arr).round(digits)

    However, this method does not always produce the same result as the current
    method, because NumPy's `round` is imprecise compared with Python's `round`
    (discussed in the NumPy documentation). For example, rounding 0.225 to two
    digits will result in 0.23 (Python) or 0.22 (NumPy).

    See Also
    --------
    util.rounder
    """
    error = 1e-7 / 10 ** digits if digits else 1e-7

    def f(num):
        near = round(num * 2, digits) / 2
        if abs(num - near) <= error:
            return round(near, digits)
        else:
            return round(num, digits)

    tmd = table.matrix_data
    tmd.data = np.vectorize(f)(tmd.data).astype('float64')
    tmd.eliminate_zeros()
    table.remove_empty(axis='observation')


def biom_add_metacol(table: biom.Table, dic, name, missing=''):
    """Add a metadata column to a table in place based on a dictionary.

    Parameters
    ----------
    table : biom.Table
        Table to add metadata column.
    dict : dict
        Metadata column (feature-to-value mapping).
    name : str
        Metadata column name.
    missing : any type, optional
        Default value if not found in dictionary.
    """
    metadata = {x: {name: dic.get(x, missing)} for x in table.ids(
        'observation')}
    table.add_metadata(metadata, axis='observation')


def clip_biom(table: biom.Table, field, sep, nested=False):
    """Clip stratified or nested feature names to a field.

    Parameters
    ----------
    table : biom.Table
        Table to clip.
    field : int
        Field index to clip at.
    sep : str
        Field separator.
    nested : bool, optional
        Whether features are nested.

    Returns
    -------
    biom.Table
        Clipped BIOM table.

    Notes
    -----
    Metadata will not be retained in the clipped table.

    See Also
    --------
    .table.clip_table
    """
    idx = field - 1

    def f1(id_, md):
        fields = id_.split(sep)
        if len(fields) >= field and fields[idx]:
            return sep.join(fields[:field]) if nested else fields[idx]

    table = table.collapse(f1, norm=False, axis='observation',
                           include_collapsed_metadata=False)
    return table.filter(lambda val, id_, md:  bool(id_), axis='observation')


def collapse_biom(table: biom.Table, mapping: dict, divide=False, field=None,
                  sep=None, nested=None):
    """Collapse a BIOM table in many-to-many mode.

    Parameters
    ----------
    table : biom.Table
        Table to collapse.
    mapping : dict of list of str
        Source-to-target(s) mapping.
    divide : bool, optional
        Whether divide per-target counts by number of targets per source.
    field : int, optional
        Index of field to be collapsed in a stratified table.
    sep : str, optional
        Field separator (must be provided if field is set).
    nested : bool, optional
        Whether features are nested.

    Returns
    -------
    biom.Table
        Collapsed BIOM table.

    Notes
    -----
    Metadata will not be retained in the collapsed table.

    See Also
    --------
    .table.collapse_table
    """
    if nested:
        f = ('{}' + sep + '{}').format
    if field:
        idx = field - 1

    # generate metadata
    metadata = {}
    for id_ in table.ids('observation'):
        feature = id_
        if field:
            fields = feature.split(sep)
            n = len(fields)
            if nested:
                fields = list(accumulate(fields, f))
            try:
                feature = fields[idx]
            except IndexError:
                continue
            if not feature:
                continue
        if feature not in mapping:
            continue
        targets = []
        for target in mapping[feature]:
            if field:
                if not nested:
                    fields[idx] = target
                    target = '|'.join(fields)
                else:
                    if idx:
                        target = fields[idx - 1] + '|' + target
                    if field < n:
                        target = target + '|' + fields[-1]
            targets.append(target)
        metadata[id_] = dict(part=targets)

    # filter table features
    table = table.filter(lambda data, id_, md: id_ in metadata,
                         axis='observation', inplace=False)

    # stop if no feature left
    if table.is_empty():
        return table

    # add mapping to table metadata
    table.add_metadata(metadata, axis='observation')

    # determine collapsing method
    kwargs = dict(norm=False, one_to_many=True, axis='observation',
                  one_to_many_mode=('divide' if divide else 'add'))

    # collapse table in many-to-many mode
    table = table.collapse(lambda _, md: zip(md['part'], md['part']), **kwargs)

    # round to integers
    if divide:
        round_biom(table)

    # clean up
    table.del_metadata(keys=['Path'])
    return table
