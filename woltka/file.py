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
from shutil import which
from subprocess import Popen, PIPE
import gzip
import bz2
import lzma


zipfmts = {'.gz':   'gzip', '.gzip':   'gzip',
           '.bz2': 'bzip2', '.bzip2': 'bzip2',
           '.xz':     'xz', '.lz':       'xz', '.lzma': 'xz'}
ziplibs = {'gzip': gzip, 'bzip2': bz2, 'xz': lzma}


def openzip(fp, mode='rt'):
    """Open a regular or compressed file by matching filename extension to
    proper library.

    Parameters
    ----------
    fp : str
        Input filepath.
    mode : str, optional
        Python file mode. Default: "rt" (read as text).

    Returns
    -------
    file handle
        Text stream ready to be read.

    Notes
    -----
    This is a simple and universal solution which uses Python's built-in
    compression modules. It supports reading and writing. However it is not as
    fast as `readzip` in reading compressed files.

    See Also
    --------
    readzip
    """
    ext = splitext(fp)[1]
    zipper = getattr(ziplibs[zipfmts[ext]], 'open') if ext in zipfmts else open
    return zipper(fp, mode)


def readzip(fp, zippers=None):
    """Open a regular or compressed file by matching filename extension to
    proper library.

    Parameters
    ----------
    fp : str
        Input filepath.
    zippers : dict of bool, optional
        Available external compression programs.

    Returns
    -------
    file handle
        Text stream ready to be read.

    Notes
    -----
    This function attempts to call external compression programs to read the
    compressed file. This is typically faster than using Python's built-in
    compression modules, because it saves the overhead in the Python wrapper,
    and utilizes additional thread(s) for decompression.

    The function checks the availability of a certain compression program at
    its first use, and stores this information in the parameter `zippers` which
    is shared across the program.

    If the external compression program is not available, or disabled (when
    `parameter` is not provided), the function will call Python's built-in
    modules instead.

    In addition to this method, a pure Python method for accelerating reading
    compression files is to use the buffered reader:

    >>> with gzip.open(fp, 'rb') as fh:
    >>>     with io.BufferedReader(fh, 1048576) as buf:
    >>>         for line in buf:
    >>>             res = line.decode()
    >>>             ...

    See Also
    --------
    openzip
    """
    # filename extension
    ext = splitext(fp)[1]

    # not a compressed file
    if ext not in zipfmts:
        return open(fp, 'r')

    fmt = zipfmts[ext]

    # external programs are disabled
    if zippers is None:
        return ziplibs[fmt].open(fp, 'rt')

    # check whether specific external program exists
    if fmt not in zippers:
        zippers[fmt] = bool(which(fmt))

    # use external program
    if zippers[fmt]:
        return Popen([fmt, '-cdfq', fp], stdout=PIPE, encoding='utf-8').stdout

    # external program does not exist
    else:
        return ziplibs[fmt].open(fp, 'rt')


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
        if ext in zipfmts:
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


def stem2rank(fp):
    """Extract rank name from stem filename.

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
    stem = path2stem(fp)

    # find pattern "a_to_b", "a-2-b" etc.
    for sep in ('-', '_'):
        try:
            left, middle, right = stem.split(sep)
        except ValueError:
            continue
        if middle in ('to', '2'):
            return right

    # find pattern "a2b"
    try:
        left, right = stem.split('2')
    except ValueError:
        pass
    else:
        return right

    # if neither is found, return entire string
    return stem


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

    Notes
    -----
    Only top-level directory is searched. Only files but not subdirectories are
    considered.
    """
    res = {}
    for fname in listdir(dir_):
        if isfile(join(dir_, fname)):
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
    with openzip(fp) as fh:
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


def read_map_uniq(fh, sep='\t'):
    """Read a mapping file, considering unique matches only (exactly two
    columns).

    Parameters
    ----------
    fh : file handle
        Mapping file.

    Yields
    ------
    tuple of (str, str)
        Key-value pair.
    """
    for line in fh:
        key, found, value = line.partition(sep)
        if found and sep not in value:
            yield key, value.rstrip()


def read_map_1st(fh, sep='\t'):
    """Read a mapping file, considering the first two columns as key and value
    while discarding remaining columns if any.

    Parameters
    ----------
    fh : file handle
        Mapping file.

    Yields
    ------
    tuple of (str, str)
        Key-value pair.
    """
    for line in fh:
        key, found, rest = line.partition(sep)
        if found:
            value, found, rest = rest.partition(sep)
            yield key, value.rstrip()


def read_map_all(fh, sep='\t'):
    """Read a mapping file, considering all columns as key (1st) and values
    (2nd-last).

    Parameters
    ----------
    fh : file handle
        Mapping file.

    Yields
    ------
    tuple of (str, list of str)
        Key-value(s) pair.
    """
    for line in fh:
        key, found, rest = line.partition(sep)
        if found:
            yield key, rest.rstrip().split(sep)


def read_map_many(fh, sep='\t'):
    """Read many-to-many relationships from a mapping file.

    Parameters
    ----------
    fh : file handle
        Mapping file.

    Returns
    -------
    dict of list
        Key-value(s) pair.

    Notes
    -----
    This is a general solution which applies to the following formats:
    1 (one-to-one):
        source1 <tab> target1
        source2 <tab> target2
        source3 <tab> target3
        ...
    2 (one-to-many):
        source1 <tab> target1
        source2 <tab> target1 <tab> target2
        source3 <tab> target2 <tab> target3 <tab> target4
        ...
    3 (many-to-many):
        source1 <tab> target1
        source1 <tab> target2
        source2 <tab> target2
        source3 <tab> target3
        source4 <tab> target3
        ...
    """
    res = {}
    for key, values in read_map_all(fh, sep):
        res.setdefault(key, []).extend(values)
    return res


def write_readmap(fh, qryque, taxque, namedic=None):
    """Write a read map to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    qryque : iterable of str
        Query sequences.
    taxque : iterable of str or dict
        Taxon(a) assigned to each query.
    namedic : dict, optional
        Taxon name dictionary.
    """
    # sort subjects by count (high-to-low) then by alphabet
    def sortkey(x): return -x[1], x[0]
    for query, taxa in zip(qryque, taxque):
        row = [query]
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
