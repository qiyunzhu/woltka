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
        elif cend >= start - 1:
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


def write_coverage(covers, outdir):
    """Write subject coverage to per sample output files.

    Parameters
    ----------
    covers : dict of dict
        Coverage data structure.
    outdir : str
        Directory of output files.
    """
    makedirs(outdir, exist_ok=True)
    for sample, cover in sorted(covers.items()):
        with open(join(outdir, f'{sample}.cov'), 'w') as fh:
            for subject, ranges in sorted(cover.items()):
                for start, end in sorted(zip(*[iter(ranges)] * 2)):
                    print(subject, start, end, sep='\t', file=fh)
