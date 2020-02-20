#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os import listdir
from os.path import join, isfile, isdir

from .util import file2stem, path2stem


def read_ids(fh):
    """Read a list of Ids from a file.

    Parameters
    ----------
    fh : file handle
        ID list file.

    Returns
    -------
    list of str
        ID list.

    Notes
    -----
    Only the first column before <tab> is considered. Lines starting with "#"
    are omitted. Empty entries are discarded.
    """
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


def match_sample_file(path_, ext=None, demux=None, samples=None):
    """Match sample IDs and alignment files.

    Parameters
    ----------
    path_ : str
        Path to a file or a directory.
    ext : str, optional
        Filename extension.
    demux : bool, optional
        Whether perform demultiplexing.
    samples : list of str, optional
        Predefined sample ID list.

    Returns
    -------
    bool, list or dict
        Whether perform demultiplexing.
        Filepaths if demultiplexing, or filepath to sample ID map if not.
    """
    errmsg = 'Provided sample IDs and actual files are inconsistent.'

    # path is an alignment file
    if isfile(path_):

        # turn on demultiplexing if not decided
        demux = demux is not False
        if demux:
            res = [path_]

        # validate with given sample Ids
        else:
            sample = path2stem(path_, ext)
            if samples and samples != [sample]:
                raise ValueError(errmsg)
            res = {path_: sample}

    # path is a directory of alignment files
    elif isdir(path_):

        # turn off demultiplexing if not decided
        demux = demux or False

        # get a map of plausible sample Ids to files
        map_ = id2file_map(path_, ext, not demux and samples)
        if len(map_) == 0:
            raise ValueError('No valid file found in directory.')
        if demux:
            res = sorted([join(path_, x) for x in map_.values()])

        # validate with given sample Ids
        else:
            if not samples:
                samples = sorted(map_.keys())
            elif len(map_) < len(samples):
                raise ValueError(errmsg)
            res = {join(path_, map_[x]): x for x in samples}

    else:
        raise ValueError(f'"{path_}" is not a valid file or directory.')

    return demux, res


def id2file_map(dir_, ext=None, ids=None):
    """Generate a map of Ids to files.

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
