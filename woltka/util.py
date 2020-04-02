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


def add_dict(dic, key, value):
    """Add a key-value pair to a dictionary while checking conflicts with
    existing ones.

    Parameters
    ----------
    dic : dict
        Dictionary to add to.
    key : hashable
        Key to be added.
    value : any
        Value to be added.

    Raises
    ------
    AssertionError
        Key already exists and its value is different.

    See Also
    --------
    update_dict

    Notes
    -----
    The behaviors of "add" and "update" (see also "update_dict") are analogous
    to Python set add() and update(), but with conflict checking.
    """
    try:
        assert dic[key] == value, f'Conflicting values found for "{key}".'
    except KeyError:
        dic[key] = value


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

    See Also
    --------
    add_dict
    sum_dict

    Notes
    -----
    If conflict-checking is not necessary, one should use Python's built-in
    method `update`.
    """
    for key, value in other.items():
        add_dict(dic, key, value)


def sum_dict(dic, other):
    """Merge one dictionary with another, adding values for common keys.

    Parameters
    ----------
    dic : dict
        Dictionary to be updated.
    other : dict
        Dictionary as source of new items.

    See Also
    --------
    update_dict
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


def feat_n_cnt(s, sep=':'):
    """Read feature and count from a string in format of "feature:count".

    Parameters
    ----------
    s : str
        Input string.
    sep : str, optional
        Separator between feature and count (default: colon).

    Returns
    -------
    tuple of (str, int)
        Pair of feature and count.
    """
    # find first occurrence of separator from right
    key, _, value = s.rpartition(sep)

    # if separator is found
    if key:

        # part after separator is number
        try:
            n = int(value)

        # otherwise, entire string is feature
        except ValueError:
            return s, 1

        # count must be positve integer
        if n > 0:
            return key, n

        # otherwise, entire string is feature
        else:
            return s, 1

    # no separator: entire string is feature
    else:
        return s, 1
