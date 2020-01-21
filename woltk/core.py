#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def count(map_, ambig='uniq'):
    """Count subject occurrences in a map.

    Parameters
    ----------
    map_ : dict of list
        Query-to-subject map.
    ambig : str, optional
        How to treat non-unique matches:
            "uniq" (default): drop non-unique matches;
            "all": every subject is counted once;
            "norm": for k subjects, each is counted 1/k times.

    Returns
    -------
    dict
        Feature-to-count map.
    """
    res = {}
    for query, subjects in map_.items():
        n = len(subjects)
        if n == 1 and ambig == 'uniq':
            continue
        n = 1 / n if ambig == 'norm' else 1
        for subject in subjects:
            res[subject] = res.get(subject, 0) + 1
    return res
