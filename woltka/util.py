#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Generic utility functions that are not specific to certain bioinformatics
operations.
"""

from os import listdir
from os.path import basename, splitext
import gzip
import bz2
import lzma


ZIPDIC = {'.gz': gzip, '.gzip': gzip,
          '.bz2': bz2, '.bzip2': bz2,
          '.xz': lzma, '.lz': lzma, '.lzma': lzma}


def readzip(fp):
    """Read a regular or compressed file by matching filename extension to
    proper library.

    Parameters
    ----------
    fp : str
        Input filepath.

    Returns
    -------
    file handle
        Text stream ready to be read.
    """
    ext = splitext(fp)[1]
    zipfunc = getattr(ZIPDIC[ext], 'open') if ext in ZIPDIC else open
    return zipfunc(fp, 'rt')


def file2stem(fname, ext=None):
    """Extract stem from filename.

    Parameters
    ----------
    fp : str
        Filepath.
    ext : str, optional
        Filename extension.

    Returns
    -------
    str
        Filename stem.

    Notes
    -----
    Common compression file extensions are recognized and stripped beforehand.

    There is no commonly accepted term for "filename with extension stripped".
    The term "stem" is from C++. In Python's `splitext` function, it is called
    "root". But "root" in this program describes tree root.
    """
    if ext is not None:
        if not fname.endswith(ext):
            raise ValueError('Filepath and filename extension do not match.')
        return fname[:-len(ext)]
    else:
        stem, ext = splitext(fname)
        if ext in ZIPDIC:
            stem, ext = splitext(stem)
        return stem


def path2stem(fp, ext=None):
    """Get filename stem from filepath.

    Parameters
    ----------
    fp : str
        Filepath.
    ext : str, optional
        Filename extension.

    Returns
    -------
    str
        Filename stem.
    """
    return file2stem(basename(fp), ext)


def update_dict(dic, other):
    """Update one dictionary with another while checking conflicting items.

    Parameters
    ----------
    dic : dict
        Dictionary to be updated.
    other : dict
        Dictionary as source of new items.

    Raises
    ------
    AssertionError
        Two dictionaries have conflicting values of the same key.

    Notes
    -----
    If conflict-checking is not necessary, one should use Python's built-in
    method `update`.
    """
    for key, value in other.items():
        try:
            assert dic[key] == value, f'Conflicting values found for "{key}".'
        except KeyError:
            dic[key] = value


def add_dict(dic, other):
    """Combine one dictionary with another, adding values for common keys.

    Parameters
    ----------
    dic : dict
        Dictionary to be updated.
    other : dict
        Dictionary as source of new items.
    """
    for key, value in other.items():
        dic[key] = dic.get(key, 0) + value


def intize(dic, zero=False):
    """Convert a dictionary of numbers to integers.

    Parameters
    ----------
    dic : dict
        Input dictionary.
    zero : bool, optional
        Whether keep zero values.
    """
    todel = []
    for key, value in dic.items():
        intval = round(value)
        if intval or zero:
            dic[key] = intval
        else:
            todel.append(key)
    for key in todel:
        del dic[key]


def delnone(dic):
    """Delete None key if any from a dictionary.

    Parameters
    ----------
    dic : dict
        Input dictionary.
    """
    try:
        del dic[None]
    except KeyError:
        pass


def allkeys(dic):
    """Get all keys in a dict of dict.

    Parameters
    ----------
    dic : dict
        Input dictionary.

    Returns
    -------
    set
        keys
    """
    return set().union(*dic.values())


def count_list(lst):
    """Count occurrences of elements of a list and return a dictionary.

    Parameters
    ----------
    lst : list
        Input list.

    Returns
    -------
    dict
        Element-to-count map.
    """
    res = {}
    for x in lst:
        res[x] = res.get(x, 0) + 1
    return res


def last_value(lst):
    """Get last value which is not None from a list.

    Parameters
    ----------
    lst : list
        Input list.

    Returns
    -------
    scalar or None
        Last element which is not None, or None if not found.
    """
    try:
        return next(x for x in reversed(lst) if x is not None)
    except StopIteration:
        pass


def read_ids(fh):
    """Read a list of IDs from a file.

    Parameters
    ----------
    fh : iterable of str
        ID list.

    Returns
    -------
    list of str
        ID list.

    Notes
    -----
    Only the first column before <tab> is considered. Lines starting with "#"
    are omitted. Empty entries are discarded.
    """
    if fh is None:
        return
    res = []
    for line in fh:
        if not line.startswith('#'):
            id_ = line.strip().split('\t', 1)[0]
            if id_:
                res.append(id_)
    if not res:
        raise ValueError('No ID is read.')
    if len(set(res)) < len(res):
        raise ValueError('Duplicate IDs found.')
    return res


def id2file_map(dir_, ext=None, ids=None):
    """Generate a map of IDs to files.

    Parameters
    ----------
    dir_ : str
        Directory containing files.
    ext : str, optional
        Filename extension.
    ids : iterable of str, optional
        Id list.

    Returns
    -------
    dict
        Id-to-file map.
    """
    res = {}
    for fname in listdir(dir_):
        try:
            id_ = file2stem(fname, ext)
        except ValueError:
            continue
        if ids and id_ not in ids:
            continue
        if id_ in res:
            raise ValueError(f'Ambiguous files for ID: "{id_}".')
        res[id_] = fname
    return res


def write_map(fh, map_, named=None):
    """Write profile to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    map_ : dict
        Map data.
    named : dict, optional
        Taxon name dictionary.
    """
    for query, taxa in map_.items():
        try:
            row = [named[query]]
        except (TypeError, KeyError):
            row = [query]
        if isinstance(taxa, dict):
            for taxon, count in taxa.items():
                row.append(f'{taxon}:{count}')
        else:
            row.append(taxa)
        print('\t'.join(row), file=fh)


def write_table(fh, data, named=None, samples=None):
    """Write profile to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    data : dict
        Profile data.
    named : dict, optional
        Taxon name dictionary.
    samples : list, optional
        Ordered sample ID list.
    """
    if samples is None:
        samples = sorted(data)
    print('#FeatureID\t{}'.format('\t'.join(samples)), file=fh)
    for key in sorted(allkeys(data)):
        # get feature name
        try:
            row = [named[key]]
        except (TypeError, KeyError):
            row = [key]
        # get feature count
        for sample in samples:
            try:
                row.append(str(data[sample][key]))
            except KeyError:
                row.append('0')
        print('\t'.join(row), file=fh)


def prep_table(profile, samples=None):
    """Convert a profile into data, index and columns, which can be further
    converted into a Pandas DataFrame or BIOM table.

    Parameters
    ----------
    profile : dict
        Input profile.

    Returns
    -------
    list of list, list, list
        Data (2D array of values).
        Index (observation Ids).
        Columns (sample Ids).
    """
    index = sorted(allkeys(profile))
    columns = samples or sorted(profile)
    data = []
    for key in index:
        row = []
        for sample in columns:
            try:
                row.append(profile[sample][key])
            except KeyError:
                row.append(0)
        data.append(row)
    return data, index, columns
