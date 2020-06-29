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


fromkeys = dict.fromkeys


def assign_none(subs, ambig=False):
    """Assign query to subjects without using a classification system.

    Parameters
    ----------
    subs : set of str
        Subjects.
    ambig : bool, optional
        Count occurrence of each subject.

    Returns
    -------
    str or dict
        Unique assignment or assignment-to-count map.
    """
    try:
        sub, = subs
        return sub
    except ValueError:
        return fromkeys(subs, 1) if ambig else None


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
    str or dict
        Unique assignment or assignment-to-count map.
    """
    try:
        sub, = subs
        return sub if subok else tree[sub]
    except ValueError:
        lca = find_lca(subs, tree)
        return None if lca == root else lca


def assign_rank(subs, rank, tree, rankdic, root=None, above=False, major=None,
                ambig=False):
    """Assign query to a fixed rank in a classification system.

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
    above : bool, optional
        Allow assignment above rank.
    major : float, optional
        Majority-rule assignment threshold.
    ambig : bool, optional
        Count occurrence of each taxon at rank.

    Returns
    -------
    str or dict
        Unique assignment or assignment-to-count map.
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
    elif ambig:
        return count_list(filter(None, taxa))
    else:
        return None


def count(matches):
    """Count occurrences of taxa in a map.

    Parameters
    ----------
    matches : dict of str or dict
        Query-to-taxon(a) map.

    Returns
    -------
    defaultdict of str: int
        Taxon-to-count map.
    """
    res = defaultdict(int)
    for taxa in matches.values():
        try:
            # unique match (scalar)
            res[taxa] += 1
        except TypeError:
            # multiple matches (dict of subject : count), to be normalized by
            # total match count
            k = 1 / sum(taxa.values())
            for taxon, n in taxa.items():
                res[taxon] += n * k
    return res


def count_strata(matches, strata):
    """Stratify taxa in a map and count occurrences.

    Parameters
    ----------
    matches : dict of str or dict
        Query-to-taxon(a) map.
    strata : dict, optional
        Read-to-feature map for stratification.

    Returns
    -------
    defaultdict of (str, str): int
        Stratified (feature, taxon)-to-count map.
    """
    res = defaultdict(int)
    for query, taxa in matches.items():
        if query in strata:
            feature = strata[query]
            if isinstance(taxa, dict):
                k = 1 / sum(taxa.values())
                for taxon, n in taxa.items():
                    taxon = (feature, taxon)
                    res[taxon] += n * k
            else:
                taxon = (feature, taxa)
                res[taxon] += 1
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
