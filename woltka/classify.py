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

from .util import count_list
from .tree import find_rank, find_lca


def assign_none(subs, uniq=False):
    """Assign query to subjects without using a classification system.

    Parameters
    ----------
    subs : set of str
        Subjects.
    uniq : bool, optional
        Assignment must be unique.

    Returns
    -------
    str or dict
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
    subs : set of str
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
    subs : set of str
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


def count_taxa(subque, taxque):
    """Count occurrences of taxa in a map.

    Parameters
    ----------
    subque : iterable of frozenset
        Subject(s) queue to manipulate.
    taxque : iterable of str or list
        Taxon(a) assigned to each query.

    Returns
    -------
    dict of dict of int
        Map of taxon to subject to count.

    See Also
    --------
    count_taxa_strat

    Notes
    -----
    Subjects (frozenset) and taxa (list) are precise matches. This is ensured
    by the iteration stability of Python sets / frozensets. See:

    .. _Source:
        https://stackoverflow.com/questions/15479928/why-is-the-order-in-
        dictionaries-and-sets-arbitrary
    """
    res = {}
    for subs, taxa in zip(subque, taxque):
        if not taxa:
            continue

        # taxa is a string (all subjects assigned to same taxon)
        try:
            res_ = res.setdefault(taxa, {})

            # only one subject (to save compute)
            try:
                sub, = subs
                res_[sub] = res_.get(sub, 0) + 1

            # multiple subjects
            except ValueError:
                k = 1 / len(subs)
                for sub in subs:
                    res_[sub] = res_.get(sub, 0) + k

        # taxa is a list (each subject corresponds to one taxon)
        except TypeError:
            k = 1 / len(list(filter(None, taxa)))
            for taxon, sub in zip(taxa, subs):
                if not taxon:
                    continue
                res_ = res.setdefault(taxon, {})
                res_[sub] = res_.get(sub, 0) + k
    return res


def count_taxa_strat(qryque, subque, taxque, strata):
    """Stratify taxa in a map and count occurrences.

    Parameters
    ----------
    qryque : iterable of str
        Query sequences.
    subque : iterable of frozenset
        Subject(s) queue to manipulate.
    taxque : iterable of str or list
        Taxon(a) assigned to each query.
    strata : dict
        Query-to-stratum map for stratification.

    Returns
    -------
    dict of dict of int
        Map of (stratum, taxon) to subject to count.

    See Also
    --------
    count_taxa
    """
    res = {}
    for query, subs, taxa in zip(qryque, subque, taxque):
        if not taxa or query not in strata:
            continue
        stratum = strata[query]
        try:
            res_ = res.setdefault((stratum, taxa), {})
            try:
                sub, = subs
                res_[sub] = res_.get(sub, 0) + 1
            except ValueError:
                k = 1 / len(subs)
                for sub in subs:
                    res_[sub] = res_.get(sub, 0) + k
        except TypeError:
            k = 1 / len(list(filter(None, taxa)))
            for taxon, sub in zip(taxa, subs):
                if not taxon:
                    continue
                res_ = res.setdefault((stratum, taxon), {})
                res_[sub] = res_.get(sub, 0) + k
    return res


def total_taxa(counts, sizes=None):
    """Sum up subject counts per taxon.

    Parameters
    ----------
    counts : dict of dict
        Counts per subject per taxon
    sizes : dict, optional
        Subject size dictionary.

    Raises
    ------
    ValueError
        Subject not found in size dictionary.

    Notes
    -----
    This function sums the subject counts per taxon into the total count of the
    taxon. If a size dictionary is provided, the count of each subject will be
    divided by its size.
    """
    if not sizes:
        for taxon, subs in counts.items():
            counts[taxon] = sum(subs.values())
    else:
        for taxon, subs in counts.items():
            try:
                counts[taxon] = sum(
                    count / sizes[sub] for sub, count in subs.items())
            except KeyError:
                raise ValueError(
                    'One or more subjects are not found in the size map.')


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
