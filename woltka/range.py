#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


"""Functions for calculating coverage of subject sequences.

The "coverage" refers to the ranges of a subject sequence (e.g., genome) that
are covered by at least one query sequence (e.g., short sequencing read).

The functions `merge_ranges` and `parse_ranges` were modified based on the
`SortedRangeList` class implemented by Daniel Hakim (@dhakim87) in the Zebra
Filter program (https://github.com/ucsd-cmi/zebra_filter).
"""

from os import makedirs
from os.path import join

from .align import iter_align


def range_mapper(fh, fmt=None, excl=None, n=1000):
    """Read an alignment file and yield maps of query to subject(s) and their
    ranges.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    fmt : str, optional
        Alignment file format.
    excl : set, optional
        Subjects to exclude.
    n : int, optional
        Number of unique queries per chunk.

    Yields
    ------
    list of str
        Query queue.
    list of dict of list of int
        Subject-to-ranges queue.

    Notes
    -----
    Same as `plain_mapper`, except that it also returns subject ranges.

    Ranges are stored as a one-dimensional, interleaved list of start1, end1,
    start2, end2, start3, end3...

    See Also
    --------
    .align.plain_mapper
    merge_ranges
    """
    it = iter_align(fh, fmt, excl, True)
    while True:
        i, done = 0, False
        qryque, subque = [None] * n, [None] * n
        for query, records in it:

            # generate a mapping of subjects to interleaved starts and ends
            ranges = {}
            for subject, _, _, start, end in records:
                ranges.setdefault(subject, []).extend((start, end))

            qryque[i] = query
            subque[i] = ranges

            i += 1
            if i == n:
                done = True
                break

        if not done:
            if i:
                yield qryque[:i], subque[:i]
            break  # pragma: no cover
        else:
            yield qryque, subque


def merge_ranges(ranges):
    """Merge short, fragmental ranges into long, continuous ones.

    Parameters
    ----------
    ranges : list of int
        Ranges to merge.

    Returns
    -------
    list of int
        Merged ranges.

    Notes
    -----
    Ranges that have overlaps will be merged into one. For example:

    >>> merge_ranges([1, 3, 2, 4, 6, 8, 7, 9])
    [1, 4, 6, 9]
    """
    res = []
    res_extend = res.extend
    cstart, cend = None, None
    for start, end in sorted(zip(*[iter(ranges)] * 2)):
        if cend is None:
            # case 1: no active range, start active range
            cstart, cend = start, end
        elif cend >= start:
            # case 2: active range continues through this range
            # extend active range
            cend = max(cend, end)
        else:
            # case 3: active range ends before this range begins
            # write new range out, then start new active range
            res_extend((cstart, cend))
            cstart, cend = start, end
    if cend is not None:
        res_extend((cstart, cend))
    return res


def parse_ranges(rmaps, covers, chunk=20000):
    """Extract range information from read maps.

    Parameters
    ----------
    rmaps : dict of dict
        Read map data structure.
    covers : dict of list
        Subject coverage data structure.
    chunk : int, optional
        Auto-compress after adding this many new ranges.

    Notes
    -----
    The `covers` data structure is a dictionary, where the key is a tuple of
    (sample, subject), and the value is a list. The first element of the list
    is an integer indicating how many new ranges have been added since last
    merge, and the remaining elements are individual ranges.

    Once the range counter reaches the threshold specified by `chunk`, an
    operation will be automatically triggered to merge small ranges into large,
    continuous ones, so as to reduce memory space and to improve downstream
    operations.

    See Also
    --------
    merge_ranges
    """
    add_cover = covers.setdefault
    for sample, (_, subque) in rmaps.items():
        for subjects in subque:
            for subject, ranges in subjects.items():
                cover = add_cover((sample, subject), [0, []])
                count = cover[0] + len(ranges)
                if count >= chunk:
                    cover[0] = 0
                    cover[1] = merge_ranges(cover[1] + ranges)
                else:
                    cover[0] = count
                    cover[1].extend(ranges)


def calc_coverage(covers):
    """Calculate coverage per sample per subject based on ranges.

    Parameters
    ----------
    covers : dict of list
        Subject coverage data structure.

    Returns
    -------
    dict of dict of list of int
        Subject coverage data.

    Notes
    -----
    This function is executed at last after all ranges have been added. It
    performs a final merge. Then it converts the dictionary with keys as
    (sample, subject) into a nested dictionary of sample -> subject -> ranges.

    See Also
    --------
    merge_ranges
    """
    res = {}
    add_sample = res.setdefault
    for (sample, subject), (_, ranges) in covers.items():
        add_sample(sample, {})[subject] = merge_ranges(ranges)
    return res


def write_coverage(covers, outdir, fmt=None):
    """Write subject coverage to per sample output files.

    Parameters
    ----------
    covers : dict of dict
        Coverage data structure.
    outdir : str
        Directory of output files.
    fmt : str, optional
        Format of output coordinates. Can be 'bed' (default) (BED-like, i.e.,
        0-based, exclusive), 'gff' (GFF-like, i.e., 1-based, inclusive), or 
        'offset(e.g., 0 or 1),i(nclusive)/e(xclusive)'.

    Notes
    -----
    BED is 0-based, exclusive (equivalent to '0,e').
    GFF is 1-based, inclusive (equivalent to '1,i').

    .. _BED format:
        https://samtools.github.io/hts-specs/BEDv1.pdf
    .. _GFF format:
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/
        gff3.md
    """

    # determine coordinate format
    begoff, endoff = 0, 0
    errmsg = f'Invalid coverage format: {fmt}.'
    if fmt is not None:
        if fmt.count(',') == 1:
            offset, code = fmt.split(',')
            try:
                begoff = int(offset)
            except ValueError:
                raise ValueError(errmsg)
            if code == 'e':
                endoff = begoff
            elif code == 'i':
                endoff = begoff - 1
            else:
                raise ValueError(errmsg)
        elif fmt.lower() == 'gff':
            begoff, endoff = 1, 0
        elif fmt.lower() != 'bed':
            raise ValueError(errmsg)

    # write coverage files
    makedirs(outdir, exist_ok=True)
    for sample, cover in sorted(covers.items()):
        with open(join(outdir, f'{sample}.cov'), 'w') as fh:
            for subject, ranges in sorted(cover.items()):
                for beg, end in sorted(zip(*[iter(ranges)] * 2)):
                    print(subject, beg + begoff, end + endoff,
                          sep='\t', file=fh)
