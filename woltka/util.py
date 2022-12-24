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


def scale_factor(s):
    """Convert scale factor string to number.

    Parameters
    ----------
    s : str
        Scale factor string.

    Returns
    -------
    int or float
        Scale factor number.

    Raises
    ------
    ValueError
        Invalid scale factor.
    """
    s = s.strip().lower()
    if s.endswith('k'):
        factor = 1000
        s = s[:-1]
    elif s.endswith('m'):
        factor = 1000000
        s = s[:-1]
    else:
        factor = 1
    try:
        return int(s) * factor
    except ValueError:
        try:
            return float(s) * factor
        except ValueError:
            raise ValueError('Invalid scale factor.')


def scale_dict(dic, factor):
    """Multiple each value of a dictionary by a scaling factor.

    Parameters
    ----------
    dic : dict
        Dictionary to be updated.
    factor : int or float
        Scaling factor.
    """
    for key in dic:
        dic[key] *= factor


def intize(num):
    """Round a floating point number to an integer.

    Parameters
    ----------
    num : float
        Number to round.

    Returns
    -------
    int
        Rounded number.

    See Also
    --------
    initize_list
    initize_dict

    Notes
    -----
    Limited by Python floating point arithmetic, odds are that when ambiguous
    assignment is on, counts will be sums of fractions and they can be rounded
    either up or down when converted to integers, depending how these numbers
    are summed. For example:

    >>> data = (1/4, 1/3, 1/4, 1/3, 1/4, 1/3, 1/3, 1/3, 1/4, 1/4, 1/4, 1/3)
    >>> total = sum(data)
    >>> total, round(total)
    (3.5, 4)

    >>> total = sum(sum(data[i:i + 4]) for i in range(0, len(data), 4))
    >>> total, round(total)
    (3.499999999999999, 3)

    To avoid this inconsistency, this function first checks if the count is
    close enough to a half number, and if yes, it will use the half number
    instead of the original number for rounding.

    Note: Multiple scientific computing languages such as Python and R's
    default round behavior is to round half numbers to the nearest even
    number. For example:

    >>> round(0.5)
    0
    >>> round(1.5)
    2
    """
    # first, round to the nearest half number
    near = round(num * 2) / 2

    # check if the difference is small enough:
    # yes - round the half number
    # no - round the original number
    if abs(num - near) <= 0.0000001:
        return round(near)
    else:
        return round(num)


def intize_list(lst):
    """Convert list elements to integers in place.

    Parameters
    ----------
    lst : list
        Input list.

    See Also
    --------
    intize

    Notes
    -----
    Use this function instead of `initize` on large lists to save some
    function-calling overhead.
    """
    for i, element in enumerate(lst):
        near = round(element * 2) / 2
        if abs(element - near) <= 0.0000001:
            lst[i] = round(near)
        else:
            lst[i] = round(element)


def intize_dict(dic, zero=False):
    """Convert dictionary values to integers in place.

    Parameters
    ----------
    dic : dict
        Input dictionary.
    zero : bool, optional
        Whether keep zero values.

    See Also
    --------
    intize

    Notes
    -----
    Use this function instead of `initize` on large dictionaries to save some
    function-calling overhead.
    """
    todel = []
    add_todel = todel.append
    for key, value in dic.items():
        near = round(value * 2) / 2
        if abs(value - near) <= 0.0000001:
            intval = round(near)
        else:
            intval = round(value)

        # keep or skip zero values
        if intval or zero:
            dic[key] = intval
        else:
            add_todel(key)
    for key in todel:
        del dic[key]


def rounder(num, digits=None):
    """Round a floating point number to given number of digits after the
    decimal point.

    Parameters
    ----------
    num : float
        Number to round.
    digits : int, optional
        Digits after the decimal point.

    Returns
    -------
    int or float
        Rounded number.

    See Also
    --------
    initize

    Notes
    -----
    Variant of `intize`, allowing rounding to custom decimal precision.
    """
    near = round(num * 2, digits) / 2
    if abs(num - near) <= (1e-7 / 10 ** digits if digits else 1e-7):
        return round(near, digits)
    else:
        return round(num, digits)


def round_list(lst, digits=None):
    """Round list elements in place.

    Parameters
    ----------
    lst : list
        Input list.
    digits : int, optional
        Digits after the decimal point.

    See Also
    --------
    rounder
    intize_list
    """
    error = 1e-7 / 10 ** digits if digits else 1e-7
    for i, element in enumerate(lst):
        near = round(element * 2, digits) / 2
        if abs(element - near) <= error:
            lst[i] = round(near, digits)
        else:
            lst[i] = round(element, digits)


def round_dict(dic, digits=None, zero=False):
    """Round dictionary values in place.

    Parameters
    ----------
    dic : dict
        Input dictionary.
    digits : int, optional
        Digits after the decimal point.
    zero : bool, optional
        Whether keep zero values.

    See Also
    --------
    rounder
    intize_dict
    """
    todel = []
    add_todel = todel.append
    error = 1e-7 / 10 ** digits if digits else 1e-7
    for key, value in dic.items():
        near = round(value * 2, digits) / 2
        if abs(value - near) <= error:
            intval = round(near, digits)
        else:
            intval = round(value, digits)
        if intval or zero:
            dic[key] = intval
        else:
            add_todel(key)
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


def feature_count(s, sep=':'):
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
    feature, _, count = s.rpartition(sep)

    # if separator is found
    if feature:

        # part after separator is number
        try:
            count = int(count)

        # otherwise, entire string is feature
        except ValueError:
            return s, 1

        # count must be positive integer
        if count > 0:
            return feature, count

        # otherwise, entire string is feature
        else:
            return s, 1

    # no separator: entire string is feature
    else:
        return s, 1
