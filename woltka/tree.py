#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for handling hierarchical classification systems.

Notes
-----
A hierarchical classification system is represented by a dict of dict which has
three attributes: "parent", "rank" and "name".

The system is based on unique identifiers (such as NCBI TaxID) instead of
descriptive names (e.g., "Escherichia coli").

In principle, the system is rank-independent (i.e., "rank" can be anything or
left empty).
"""

from woltka.util import readzip


# standard taxonomic ranks
rankorder = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus',
             'species', 'strain')
rank2code = {'kingdom': 'k', 'phylum': 'p', 'class':   'c', 'order':  'o',
             'family':  'f', 'genus':  'g', 'species': 's', 'strain': 't'}
code2rank = {v: k for k, v in rank2code.items()}

# unclassified taxon names
notax = {'', '0', 'unclassified', 'unassigned'}


def build_tree(names_fp: str = None,
               nodes_fp: str = None,
               newick_fp: str = None,
               levels_fp: str = None,
               lineage_fp: str = None,
               group_fps: list = None) -> (dict, dict, dict):
    """Construct a tree to represent the hierarchical classification system.

    Parameters
    ----------
    names_fp : str
        Taxonomic names file.
    nodes_fp : str
        Taxonomic nodes file.
    newick_fp : str
        Newick tree file.
    levels_fp : str
        Taxonomic level table file.
    lineage_fp : str
        Lineage strings file.
    group_fps : list of str
        Grouping files in order.

    Returns
    -------
    dict, dict, dict
        Taxonomic tree.
        Taxon rank map.
        Taxon name map.
    """
    tree, ranks, names = {}, {}, {}
    if names_fp:
        with readzip(names_fp) as f:
            names = read_names(f)
    if nodes_fp:
        with readzip(nodes_fp) as f:
            tree, ranks = read_nodes(f)
        names = {k: v for k, v in names.items() if k in tree}
        return tree, ranks, names
    if newick_fp:
        with readzip(newick_fp) as f:
            tree = read_newick(f)
        names = {k: v for k, v in names.items() if k in tree}
        return tree, ranks, names
    if levels_fp:
        with readzip(levels_fp) as f:
            tree = read_levels(f)
        names = {k: v for k, v in names.items() if k in tree}
        return tree, ranks, names


def read_names(fh):
    """Read taxonomic names from a file.

    Parameters
    ----------
    fh : file handle
        Taxon name mapping file.

    Returns
    -------
    dict
        Taxon name map.

    Notes
    -----
    Can be NCBI-style names.dmp or a plain map of Id to name.
    """
    names = {}
    for line in fh:
        x = line.rstrip('\r\n').replace('\t|', '').split('\t')
        if len(x) < 4 or x[3] == 'scientific name':
            names[x[0]] = x[1]
    return names


def read_nodes(fh):
    """Read taxonomic nodes from a file.

    Parameters
    ----------
    fh : file handle
        Taxonomic nodes file.

    Returns
    -------
    dict, dict
        Taxonomy tree, taxonomic rank map.

    Notes
    -----
    Input file can be NCBI-style nodes.dmp or a plain map of Id to parent and
    (optional) rank.
    """
    tree, ranks = {}, {}
    for line in fh:
        x = line.rstrip('\r\n').replace('\t|', '').split('\t')
        tree[x[0]] = x[1]
        try:
            ranks[x[0]] = x[2]
        except IndexError:
            pass
    return tree, ranks


def read_newick(fh):
    return


def read_levels(fh):
    return


def read_lineage(fh):
    """Read taxonomic lineage strings from a file and build a taxonomy tree.

    Parameters
    ----------
    fh : file handle
        Taxonomic lineage mapping file.

    Returns
    -------
    dict
        Taxonomy tree.

    Raises
    ------
    ValueError
        Conflict of taxon-parent relationship is found.

    Notes
    -----
    A lineage string consists of one or multiple semicolon-delimited taxa.
    Spaces leading or trailing a taxon are stripped. Spaces within a taxon
    are tolerable. Each taxon must be unique.
    """
    tree = {}
    for line in fh:
        id_, lineage = line.rstrip('\r\n').split('\t')
        parent = None
        for taxon in lineage.split(';'):
            # strip white spaces
            taxon = taxon.strip()
            # skip empty taxon
            # if taxon.lower() in notax or taxon[1:] == '__':
            #     continue
            this = f'{parent};{taxon}' if parent else taxon
            if this not in tree:
                tree[this] = parent
            elif tree[this] != parent:
                raise ValueError(f'Conflict in taxon "{this}".')
            parent = this
    return tree
