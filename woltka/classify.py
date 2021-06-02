#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for classifying query sequences by assigning their subjects to a
hierarchical classification system.
"""

from operator import itemgetter
from collections import defaultdict

from .util import count_list
from .tree import find_rank, find_lca


"""
Assigner: assign a query sequence to one or more taxa based on its subjects and
a classification system, if applicable.

- assign_none
- assign_free
- assign_rank
"""


def assign_none(subs, uniq=False):
    """Assign query to subjects without using a classification system.

    Parameters
    ----------
    subs : tuple of str
        Subjects.
    uniq : bool, optional
        Assignment must be unique.

    Returns
    -------
    str or list or None
        Unique subject or list of subjects.
    """
    try:
        sub, = subs
        return sub
    except ValueError:
        return None if uniq else list(subs)


def assign_free(subs, tree, root=None, subok=False):
    """Assign query based on a classification system in a rank-free manner.

    Parameters
    ----------
    subs : tuple of str
        Subjects.
    tree : dict
        Hierarchical classification system.
    root : str, optional
        Root identifier.
    subok : bool, optional
        Allow assignment to subjects.

    Returns
    -------
    str or None
        Unique assignment.
    """
    try:
        sub, = subs
        return sub if subok else (tree[sub] if sub in tree else None)
    except ValueError:
        lca = find_lca(subs, tree)
        return None if lca == root else lca


def assign_rank(subs, rank, tree, rankdic, root=None, major=None, above=False,
                uniq=False):
    """Assign query to a given rank in a classification system.

    Parameters
    ----------
    subs : tuple of str
        Subjects.
    rank : str
        Target rank.
    tree : dict
        Hierarchical classification system.
    rankdic : dict
        Rank dictionary.
    root : str, optional
        Root identifier.
    major : float, optional
        Majority-rule assignment threshold.
    above : bool, optional
        Allow assignment above rank.
    uniq : bool, optional
        Assignment must be unique.

    Returns
    -------
    str or list or None
        Unique assignment or list of assignments.

    TODO
    ----
    Combine major and above.
    """
    taxa = [find_rank(x, rank, tree, rankdic) for x in subs]
    tset = set(taxa)
    if len(tset) == 1:
        return taxa[0]
    elif major:
        return majority(taxa, major)
    elif above:
        if None in tset:
            return None
        lca = find_lca(tset, tree)
        return None if lca == root else lca
    elif uniq:
        return None
    else:
        return taxa


"""
Counter: convert a read-to-taxa map into a taxon-to-count map.

Four counter functions are implemented to tackle alternative parameter settings
(normalization by subject size and stratification by query), while maintaining
performance by saving conditional statements.

- counter
- counter_size
- counter_strat
- counter_size_strat
"""


def counter(taxque):
    """Count occurrences of taxa in a read-to-taxa map.

    Parameters
    ----------
    taxque : iterable of str or list

    Returns
    -------
    dict of int
        Map of taxon to count.
    """
    res = defaultdict(int)
    for taxa in taxque:
        if not taxa:
            continue

        # taxa is string: unique match
        try:
            res[taxa] += 1

        # taxa is list: multiple matches, to be divided by match count
        except TypeError:
            taxa = list(filter(None, taxa))
            k = 1 / len(taxa)
            for taxon in taxa:
                res[taxon] += k
    return res


def counter_size(subque, taxque, sizes):
    """Count occurrences of taxa in a read map while normalizing by subject
    size.

    Parameters
    ----------
    subque : iterable of tuple
        Subject(s) mapped to each query.
    taxque : iterable of str or list
        Taxon(a) assigned to each query.
    sizes : dict
        Subject size dictionary.

    Returns
    -------
    dict of float
        Map of taxon to normalized count.

    Raises
    ------
    KeyError
        Subject not found in size dictionary.

    See Also
    --------
    counter
    """
    res = defaultdict(int)
    for subs, taxa in zip(subque, taxque):
        if not taxa:
            continue
        try:
            res[taxa] += sum(sizes[x] for x in subs) / len(subs)
        except TypeError:
            k = 1 / len(list(filter(None, taxa)))
            for taxon, sub in zip(taxa, subs):
                if not taxon:
                    continue
                res[taxon] += sizes[sub] * k
    return res


def counter_strat(qryque, taxque, strata):
    """Stratify taxa in a read map and count occurrences.

    Parameters
    ----------
    qryque : iterable of str
        Query sequences.
    taxque : iterable of str or list
        Taxon(a) assigned to each query.
    strata : dict
        Query-to-stratum map.

    Returns
    -------
    dict of int
        Map of (stratum, taxon) to count.

    See Also
    --------
    counter
    """
    res = defaultdict(int)
    for query, taxa in zip(qryque, taxque):
        if not taxa or query not in strata:
            continue
        stratum = strata[query]
        try:
            res[(stratum, taxa)] += 1
        except TypeError:
            taxa = list(filter(None, taxa))
            k = 1 / len(taxa)
            for taxon in taxa:
                res[(stratum, taxon)] += k
    return res


def counter_size_strat(qryque, subque, taxque, sizes, strata):
    """Stratify taxa in a read map and count occurrences while normalizing by
    subject size.

    Parameters
    ----------
    qryque : iterable of str
        Query sequences.
    subque : iterable of tuple
        Subject(s) mapped to each query.
    taxque : iterable of str or list
        Taxon(a) assigned to each query.
    sizes : dict, optional
        Subject size dictionary.
    strata : dict
        Query-to-stratum map.

    Returns
    -------
    dict of float
        Map of (stratum, taxon) to normalized count.

    Raises
    ------
    KeyError
        Subject not found in size dictionary.

    See Also
    --------
    counter_size
    counter_strat
    """
    res = defaultdict(int)
    for query, subs, taxa in zip(qryque, subque, taxque):
        if not taxa or query not in strata:
            continue
        stratum = strata[query]
        try:
            res[(stratum, taxa)] += sum(sizes[x] for x in subs) / len(subs)
        except TypeError:
            k = 1 / len(list(filter(None, taxa)))
            for taxon, sub in zip(taxa, subs):
                if not taxon:
                    continue
                res[(stratum, taxon)] += sizes[sub] * k
    return res


def majority(taxa, th=0.8):
    """Select taxon from list by majority rule.

    Parameters
    ----------
    taxa : list of str
        Input taxon list.
    th : float, optional
        Threshold of majority, range = (0.5, 1.0].

    Returns
    -------
    str or None
        Selected taxon.
    """
    for taxon, n in sorted(count_list(taxa).items(), key=itemgetter(1),
                           reverse=True):
        return taxon if n >= len(taxa) * th else None
