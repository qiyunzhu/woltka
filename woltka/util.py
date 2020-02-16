#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

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
        input filepath

    Returns
    -------
    file handle
        text stream ready to be read
    """
    ext = splitext(fp)[1]
    zipfunc = getattr(ZIPDIC[ext], 'open') if ext in ZIPDIC else open
    return zipfunc(fp, 'rt')


def path2stem(fp):
    """Get filename stem from filepath.

    Parameters
    ----------
    fp : str
        filepath

    Returns
    -------
    str
        filename stem
    """
    stem, ext = splitext(basename(fp))
    if ext in ZIPDIC:
        stem, ext = splitext(stem)
    return stem


def id2file_map(dir_, ext=None, ids=None):
    """Generate a map of Ids to files.

    Parameters
    ----------
    dir_ : str
        directory to search for matching files
    ext : str, optional
        filename extension
    ids : iterable of str, optional
        Id list

    Returns
    -------
    dict
        Id to file map
    """
    res = {}
    if ext is not None:
        if not ext.startswith('.'):
            ext = '.' + ext
        n = len(ext)
    for fname in listdir(dir_):
        id_ = None
        if ext is not None:
            if fname.endswith(ext):
                id_ = fname[:-n]
        else:
            id_, ext_ = splitext(fname)
            if ext_ in ZIPDIC:
                id_ = splitext(id_)[0]
        if id_ is None:
            continue
        if ids is not None and id_ not in ids:
            continue
        if id_ in res:
            raise ValueError('Ambiguous files for Id: {}.'.format(id_))
        res[id_] = fname
    return res


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


def integerize(dic):
    """Convert a map of numbers to integers.

    Parameters
    ----------
    dic : dict
        input dictionary

    Returns
    -------
    dict
        converted dictionary
    """
    return {k: int(v) for k, v in dic.items() if int(v)}


def allkeys(dic):
    """Get all keys in a dict of dict.

    Parameters
    ----------
    dic : dict
        input dictionary

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
