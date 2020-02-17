#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Core functionalities.
"""

from woltka.util import count_list
from woltka.tree import find_lca, find_rank


def count(matches):
    """Count occurrences of taxa in a map.

    Parameters
    ----------
    matches : dict of str or dict
        Query-to-taxon(a) map.

    Returns
    -------
    dict
        Taxon-to-count map.
    """
    res = {}
    for taxa in matches.values():
        try:
            # unique match (scalar)
            res[taxa] = res.get(taxa, 0) + 1
        except TypeError:
            # multiple matches (dict of subject : count), to be normalized by
            # total match count
            k = 1 / sum(taxa.values())
            for taxon, n in taxa.items():
                res[taxon] = res.get(taxon, 0) + n * k
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
    for taxon, n in sorted(count_list(taxa).items(), key=lambda x: x[1],
                           reverse=True):
        return taxon if n >= len(taxa) * th else None


def assign(subs, rank=None, tree=None, rankd=None, root=None,
           above=False, multi=False, major=None, subok=False):
    """Assign a query sequence to a classification unit based on its subjects.

    Parameters
    ----------
    subs : set of str
        Subject(s) of a query sequence.
    rank : str, optional
        Ranks to assign to, or "free" for rank-free LCA assignment.
    tree : dict, optional
        Hierarchical classification system.
    rankd : dict, optional
        Rank dictionary.
    root : str, optional
        Root identifier.
    above : bool, optional
        Assignment above given rank is acceptable (for fixed ranks).
    major : float, optional
        Majority-rule assignment threshold (available only with a fixed rank
        and not above or multi).
    multi : bool, optional
        Count occurrence of each possible assignment instead of targeting one
        assignment (available only with a fixed rank and not above).
    subok : bool, optional
        Allow assignment to subject(s) itself instead of higher classification
        units.

    Returns
    -------
    str or dict
        Unique assignment or assignment-to-count map.
    """
    # no classification, just subject(s) itself
    if rank is None or rank == 'none' or tree is None:
        if len(subs) == 1:
            return max(subs)
        elif multi:
            return count_list(subs)
        else:
            return None

    # free rank classification: find LCA
    elif rank == 'free':
        if len(subs) == 1:
            if subok:
                return tree[max(subs)]
            else:
                return max(subs)
        else:
            lca = find_lca(subs, tree)
            return None if lca == root else lca

    # fixed rank classification
    else:
        taxa = [find_rank(x, rank, tree, rankd) for x in subs]
        if len(set(taxa)) == 1:
            return taxa[0]
        elif major:
            return majority(taxa, major)
        elif above:
            lca = find_lca(set(taxa), tree)
            return None if lca == root else lca
        elif multi:
            return count_list(taxa)
        else:
            return None
