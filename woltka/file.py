#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for handling input and output files.
"""

from os import listdir
from os.path import basename, dirname, splitext, isfile, join
import gzip
import bz2
import lzma

from .util import allkeys
from .tree import get_lineage_gg


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

    Raises
    ------
    ValueError
        Multiple files have the same stem filename.
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
            fp = join(fdir, file_)
            if isfile(fp):
                res.append((id_, fp))
                continue

            # if previous lines appear to be valid files, raise error
            if res:
                raise ValueError(f'Alignment file "{file_}" does not exist.')

            # otherwise (i.e., this is the first line), return
            return
    return res or None


def write_readmap(fh, rmap, namedic=None):
    """Write a read map to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    rmap : dict
        Read-to-taxon(a) map.
    namedic : dict, optional
        Taxon name dictionary.
    """
    # sort subjects by count (high-to-low) then by alphabet
    def sortkey(x): return -x[1], x[0]
    for read, taxa in rmap.items():
        row = [read]
        if isinstance(taxa, dict):
            for taxon, count in sorted(taxa.items(), key=sortkey):
                if namedic and taxon in namedic:
                    taxon = namedic[taxon]
                row.append(taxon + ':' + str(count))
        elif namedic and taxa in namedic:
            row.append(namedic[taxa])
        else:
            row.append(taxa)
        print('\t'.join(row), file=fh)


def write_table(fh, data, samples=None, tree=None, rankdic=None, namedic=None,
                name_as_id=False):
    """Write a profile to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    data : dict
        Profile data.
    samples : list, optional
        Ordered sample ID list.
    tree : dict, optional
        Taxonomic tree, to inform "Lineage" column.
    rankdic : dict, optional
        Rank dictionary, to inform "Rank" column.
    namedic : dict, optional
        Taxon name dictionary, to inform "Name" column.
    name_as_id : bool, optional
        Replace feature IDs with names. It applies to row headers and "Lineage"
        column, and removes "Name" column.

    Notes
    -----
    The output table will have columns as samples and rows as features.
    Optionally, three metadata columns, "Name", "Rank" and "Lineage" will be
    appended to the right of the table.
    """
    samples = samples or sorted(data)

    # table header
    header = ['#FeatureID'] + samples
    if namedic and not name_as_id:
        header.append('Name')
    if rankdic:
        header.append('Rank')
    if tree:
        header.append('Lineage')
    print('\t'.join(header), file=fh)

    # table body
    for key in sorted(allkeys(data)):
        # get feature name
        name = namedic[key] if namedic and key in namedic else None
        # fill row header (feature Id or name)
        row = [namedic[key]] if name_as_id and name else [key]
        # fill cell values (feature counts)
        for sample in samples:
            row.append(str(data[sample][key]) if key in data[sample] else '0')
        # fill name column
        if namedic and not name_as_id:
            row.append(name or '')
        # fill rank column
        if rankdic:
            row.append(rankdic[key] if key in rankdic else '')
        # fill lineage column
        if tree:
            row.append(get_lineage_gg(
                key, tree, namedic if name_as_id else None))
        # print row
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

    Notes
    -----
    This function is currently not in use.
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
