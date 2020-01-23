#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os import listdir
from os.path import splitext
import gzip
import bz2
import lzma


zipdict = {'.gz': gzip, '.gzip': gzip,
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
    zipfunc = getattr(zipdict[ext], 'open') if ext in zipdict else open
    return zipfunc(fp, 'rt')


def id2file_map(dir_, ext=None, ids=None):
    """Generate a map of Ids to files.

    Parameters
    ----------
    dir_ : str
        directory to search for matching files
    ext : str, optional
        filename extension
    ids : iterable, optional
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
            if ext_ in zipdict:
                id_ = splitext(id_)[0]
        if id_ is None:
            continue
        if ids is not None and id_ not in ids:
            continue
        if id_ in res:
            raise ValueError('Ambiguous files for Id: {}.'.format(id_))
        res[id_] = fname
    return res


def integerize(dic):
    """Convert a map of numbers to integers."""
    return {k: int(v) for k, v in dic.items() if int(v)}


def allkeys(dic):
    """Get all keys in a dict of dict."""
    return set().union(*dic.values())
