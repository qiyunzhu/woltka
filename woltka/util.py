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
from os.path import basename, dirname, splitext, isfile, join
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
        Keys.
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
            id_ = line.strip().partition('\t')[0]
            if id_:
                res.append(id_)
    if not res:
        raise ValueError('No ID is read.')
    if len(set(res)) < len(res):
        raise ValueError('Duplicate IDs found.')
    return res


def id2file_from_dir(dir_, ext=None, ids=None):
    """Generate an ID-to-file map from a directory.

    Parameters
    ----------
    dir_ : str
        Directory containing files.
    ext : str, optional
        Filename extension.
    ids : iterable of str, optional
        ID list.

    Returns
    -------
    dict
        ID-to-file map.
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


def id2file_from_map(fp):
    """Read ID-to-file mapping from a file.

    Parameters
    ----------
    fp : str
        Path to an ID-to-file mapping file (format: ID <tab> filepath).

    Returns
    -------
    list of tuple, or None
        Ordered ID-to-file map, or None if not a valid mapping file.

    Raises
    ------
    ValueError
        If any of 2nd to nth filepaths does not exist.

    Notes
    -----
    This function checks every line of a file, except for empty lines and lines
    starting with "#", to see if it suffices the format of sample ID <tab>
    filepath.

    The filepath must point to an existing file. The function first checks the
    filepath as absolute, or relative to the current directory. If not found,
    it then searches the same directory as the mapping file. Therefore, instead
    of specifying full paths, one could also provide just filenames as long as
    the mapping file and the alignment files are in the same directory.

    If the filepath in the first line does not exist, the function will return
    None, assuming this is not a mapping file. However if the first filepath
    exists whereas any of the filepaths below does not exist, an error will be
    raised, reminding user of the potentially incorrect filepath.
    """
    res = []
    fdir = dirname(fp)
    with readzip(fp) as fh:
        for line in fh:
            line = line.rstrip()

            # skip empty or commented lines
            if not line or line.startswith('#'):
                continue

            # a valid map must have exactly two columns per line
            try:
                id_, file_ = line.split('\t')
            except ValueError:
                return

            # check full filepath
            if isfile(file_):
                res.append((id_, file_))
                continue

            # search same directory
            path_ = join(fdir, file_)
            if isfile(path_):
                res.append((id_, path_))
                continue

            # if previous lines appear to be valid files, raise error
            if res:
                raise ValueError(f'Alignment file "{file_}" does not exist.')

            # otherwise (i.e., this is the first line), return
            return
    return res or None


def write_map(fh, map_, namedic=None):
    """Write a read map to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    map_ : dict
        Map data.
    namedic : dict, optional
        Taxon name dictionary.
    """
    for query, taxa in map_.items():
        try:
            row = [namedic[query]]
        except (TypeError, KeyError):
            row = [query]
        if isinstance(taxa, dict):
            for taxon, count in taxa.items():
                row.append(f'{taxon}:{count}')
        else:
            row.append(taxa)
        try:
            print('\t'.join(row), file=fh)
        except TypeError:
            raise ValueError(row)


def write_table(fh, data, namedic=None, samples=None):
    """Write a profile to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    data : dict
        Profile data.
    namedic : dict, optional
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
            row = [namedic[key]]
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
    list of list
        Data (2D array of values).
    list
        Index (observation Ids).
    list
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
