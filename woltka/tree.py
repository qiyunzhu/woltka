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
A hierarchical classification system is represented by a tree, which is a dict
of taxon ID to parent taxon ID. Each taxon ID is a unique identifier (such as
NCBI TaxID), although a descriptive name (e.g., "Escherichia coli") may work as
long as it is unique. The root of the tree is indicated by ID == parent ID. The
terminal leaves of the tree are IDs of subjects to which input data are mapped.
This simple data structure is sufficient for the purpose of this program.

In addition, the name and rank assignments of taxa are recorded in two separate
dicts, if available. In principle, the system is rank-independent (i.e., "rank"
can be anything or left empty).
"""

import re

from .util import readzip, path2stem, update_dict, last_value


# standard taxonomic ranks
rankorder = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus',
             'species', 'strain')
rank2code = {'kingdom': 'k', 'phylum': 'p', 'class':   'c', 'order':  'o',
             'family':  'f', 'genus':  'g', 'species': 's', 'strain': 't'}
code2rank = {v: k for k, v in rank2code.items()}

# compatibilities for top-level ranks (domain, superkingdom, kingdom)
rank2code['domain'] = 'k'
rank2code['superkingdom'] = 'k'
code2rank['d'] = 'kingdom'

# unclassified taxon names
notax = {'', '0', 'unclassified', 'unassigned'}


def build_tree(map_fps, names_fp, nodes_fp, newick_fp, ranktb_fp, lineage_fp):
    """Construct a tree to represent the hierarchical classification system.

    Parameters
    ----------
    map_fps : list of str
        Mapping files.
    names_fp : str
        Taxonomic names file.
    nodes_fp : str
        Taxonomic nodes file.
    newick_fp : str
        Newick tree file.
    ranktb_fp : str
        Rank table file.
    lineage_fp : str
        Lineage strings file.

    Returns
    -------
    dict, dict, dict, str
        Taxonomic tree.
        Rank dictionary.
        Name dictionary.
        Root identifier.
    """
    tree, rankd, named = {}, {}, {}

    # plain maps
    for fp in map_fps:
        rank = path2stem(fp)  # filename stem as rank
        with readzip(fp) as f:
            map_ = read_map(f)
        update_dict(tree, map_)
        update_dict(rankd, {k: rank for k in set(map_.values())})

    # taxdump-style names and nodes
    if names_fp:
        with readzip(names_fp) as f:
            update_dict(named, read_names(f))
    if nodes_fp:
        with readzip(nodes_fp) as f:
            tree_, rankd_ = read_nodes(f)
        update_dict(tree, tree_)
        update_dict(rankd, rankd_)

    # Newick-format tree
    if newick_fp:
        with readzip(newick_fp) as f:
            update_dict(tree, read_newick(f))

    # rank table file
    if ranktb_fp:
        with readzip(ranktb_fp) as f:
            update_dict(tree, read_ranktb(f))

    # lineage strings file
    if lineage_fp:
        with readzip(lineage_fp) as f:
            tree_, rankd_ = read_lineage(f)
        update_dict(tree, tree_)
        update_dict(rankd, rankd_)

    # fill root
    root = fill_root(tree)

    return tree, rankd, named, root


def read_map(fh):
    """Read simple low-to-high map from file.

    Parameters
    ----------
    fh : file handle
        Simple mapping file.

    Returns
    -------
    dict
        Mapping (mutually identical to a taxonomy tree).

    Notes
    -----
    Only first two columns are considered.

    TODO
    ----
    Support for one-to-multiple mappings.
    """
    try:
        return dict(x.split('\t', 2)[:2] for x in fh.read().splitlines())
    except ValueError:
        raise ValueError('Invalid mapping file format.')


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
        x = line.rstrip().replace('\t|', '').split('\t')
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
    tree, rankd = {}, {}
    for line in fh:
        x = line.rstrip().replace('\t|', '').split('\t')
        tree[x[0]] = x[1]
        try:
            rankd[x[0]] = x[2]
        except IndexError:
            pass
    return tree, rankd


def read_newick(fh):
    """Read tree from a Newick-format file.

    Parameters
    ----------
    fh : file handle
        Newick tree file.

    Returns
    -------
    dict
        Taxonomy tree.

    Raises
    ------
    ValueError
        Missing internal node Id.
        Found non-unique node Id.

    Notes
    -----
    This simple Newick parser only captures one piece of information: child-to-
    parent mapping. All others (such as branch lengths) are omitted. It assumes
    that the Newick format is correct, and only raises at two scenarios where
    the tree does not satisfy the need of this program.
    """
    # get clean Newick string
    nwk = ''.join(x.strip() for x in fh).rstrip(';')
    res = {}

    # Newick node ("(a,b)")
    pnod = re.compile(r'\([^()]+\)')

    # Newick node end ("," or ")")
    pend = re.compile(r'[,)]')

    # Newick node label ("Id:length")
    def _get_id(label):
        return label.split(':', 1)[0].strip('"\'')

    # recursively search for nodes from inside to outside
    while True:
        m = pnod.search(nwk)
        if m is None:
            break

        # get parent Id of current node
        tail = nwk[m.end(0):]
        parent = _get_id(pend.split(tail, 1)[0])
        if parent == '':
            raise ValueError('Missing internal node Id.')

        # get child Ids of current node
        for child in [_get_id(x) for x in m.group(0)[1:-1].split(',')]:
            if child in res:
                raise ValueError(f'Found non-unique node Id: "{child}".')
            res[child] = parent

        # remove node from search string
        nwk = nwk[:m.start(0)] + tail

    # outest-most parent must be root
    res[parent] = parent
    return res


def read_ranktb(fh):
    """Read taxonomic information from a rank table.

    Parameters
    ----------
    fh : file handle
        Taxonomic ranks file.

    Returns
    -------
    dict, dict
        Taxonomy tree, rank dictionary.

    Notes
    -----
    Taxonomic units must be unique.

    TODO
    ----
    Add support for same taxon name but different ranks. This may be useful.
    For example, "Actinobacteria" is both a phylum and a class.
    """
    tree, rankd = {}, {}
    ranks = None
    for line in fh:
        row = line.rstrip().split('\t')

        # get rank names from header
        if ranks is None:
            ranks = row[1:]
            continue

        # get lineage (taxa from high to low)
        lineage = [None if x in notax else x for x in row[1:]]

        # map entry to lowest taxon
        tree[row[0]] = last_value(lineage)

        for i, taxon in enumerate(lineage):
            if taxon is None:
                continue

            # map lower taxon to higher taxon
            rank, parent = ranks[i], last_value(lineage[:i])

            # check existing relationship
            try:
                if tree[taxon] != parent or rankd[taxon] != rank:
                    raise ValueError(f'Conflict at taxon "{taxon}".')

            # add current taxon to dictionary
            except KeyError:
                tree[taxon], rankd[taxon] = parent, rank

    return tree, rankd


def read_lineage(fh):
    """Read taxonomic information from a lineage mapping file.

    Parameters
    ----------
    fh : file handle
        Taxonomic lineage mapping file.

    Returns
    -------
    dict, dict
        Taxonomy tree, rank dictionary.

    Raises
    ------
    ValueError
        Conflict of taxon-parent relationship is found.

    Notes
    -----
    The lineage mapping file format is a.k.a. Greengenes-style.

    A lineage string consists of one or multiple semicolon-delimited taxa.
    Spaces leading or trailing a taxon are stripped. Spaces within a taxon
    are tolerable.

    Each taxon is represented by the entire ancestral lineage before it.
    Therefore, the taxon itself does not have to be unique, but the ancestral
    lineage does.

    Empty taxon is allowed, in consistency with QIIME convention (e.g.,
    "k__Bacteria;p__" is considered a valid phylum).
    """
    p = re.compile(r'([a-z])__.*')
    tree, rankd = {}, {}
    for line in fh:
        if line.startswith('#'):
            continue
        id_, lineage = line.rstrip().split('\t')
        parent = None
        for taxon in lineage.split(';'):
            taxon = taxon.strip()

            # skip empty taxon (currently disabled)
            # if taxon.lower() in notax or taxon[1:] == '__':
            #     continue

            # append entire ancestral lineage
            this = f'{parent};{taxon}' if parent else taxon

            # check existing relationship
            try:
                if tree[this] != parent:
                    raise ValueError(f'Conflict at taxon "{this}".')

            # add current taxon to dictionary
            except KeyError:
                tree[this] = parent

                # get rank from prefix (e.g., "k__")
                try:
                    rankd[this] = code2rank[p.match(taxon).group(1)]
                except (AttributeError, KeyError):
                    pass

            parent = this

        # add entry to whole lineage
        tree[id_] = parent

    return tree, rankd


def fill_root(tree):
    """Add a root node to a tree if absent.

    Parameters
    ----------
    tree : dict
        Taxonomy tree.

    Returns
    -------
    str
        Root node identifier

    Notes
    -----
    A root is defined as having parent as itself, a behavior derived from the
    NCBI convention. Only root must be present in a tree, so that all taxa can
    be traced back to the same root.

    In custom trees, there may or may not be a clearly defined root. This
    function aims as defining a root for any given tree. Specifically, if
    there is one root node (highest hierarchy), this node will be "sealed"
    (making parent as itself). If there are multiple crown nodes ("crown"
    describes the root of a clade instead of the entire tree), a new root
    node will be generated and serve as the parent of all crown nodes.
    """
    crown, toadd, tested = [], set(), set()
    for taxon in tree:
        this = taxon
        while True:
            if this in tested:
                break
            tested.add(this)

            # parent is missing
            try:
                parent = tree[this]
            except KeyError:
                crown.append(this)
                toadd.add(this)
                break

            # parent is None
            if parent is None:
                crown.append(this)
                break

            # parent is itself
            if parent == this:
                crown.append(this)
                break
            this = parent

    # fill non-existent root or crown nodes
    for node in toadd:
        tree[node] = None

    # this happens only when tree is empty
    if len(crown) == 0:
        return None

    # there is only one crown node
    elif len(crown) == 1:

        # make the parent of root iself
        root = crown[0]
        tree[root] = root
        return root

    # there are more than one crown node
    else:

        # in NCBI convention, root should have identifier "1"
        i = 1

        # in case "1" is already in tree, find an unused integer
        while True:
            if str(i) not in tree:
                break
            i += 1
        root = str(i)
        tree[root] = root

        # coalesce all crown nodes to root
        for x in crown:
            tree[x] = root
        return root


def get_lineage(taxon, tree):
    """Get lineage of given taxon in taxonomy tree.

    Parameters
    ----------
    taxon : str
        Query taxon.
    tree : dict
        Taxonomy tree.

    Returns
    -------
    list of str
        Lineage from root to query taxon.
    """
    # if taxon is not in tree, return None
    try:
        parent = tree[taxon]
    except KeyError:
        return None

    # move up classification hierarchy till root
    lineage = [taxon]
    this = taxon
    while True:
        lineage.append(parent)
        this = parent
        parent = tree[this]
        if parent == this:
            break

    # reverse lineage so that it's high-to-low
    return list(reversed(lineage))


def find_rank(taxon, rank, tree, rankd):
    """Find ancestor at a given rank by stepping up the classification
    hierarchy.

    Parameters
    ----------
    taxon : str
        Query taxon.
    rank : str
        Target rank.
    tree : dict
        Taxonomy tree.
    rankd : dict
        Taxon-to-rank map.

    Returns
    -------
    str
        Ancestral taxon if found.
    """
    # if taxon is not in tree, return None
    try:
        parent = tree[taxon]
    except KeyError:
        return None

    # move up hierarchy until reaching given rank
    this = taxon
    while True:

        # check rank of current taxon
        try:
            if rankd[this] == rank:
                return this
        except KeyError:
            pass

        # stop when reaching root
        if parent == this:
            break

        # move up to parent
        this = parent
        parent = tree[this]


def find_lca(taxa, tree):
    """Find lowest common ancestor (LCA) of given taxa.

    Parameters
    ----------
    taxa : iterable of str
        Query taxa.
    tree : dict
        Taxonomy tree.

    Returns
    -------
    str or None
        LCA of taxa, or None if any taxon is not in tree.

    TODO
    ----
    Combine LCA and majority rule, which is not trivial and requires careful
    reasoning and algorithm design.
    """
    taxa_ = list(taxa)

    # get lineage of first taxon
    lineage = get_lineage(taxa_[0], tree)
    if lineage is None:
        return

    # compare with remaining taxa
    for taxon in taxa_[1:]:
        try:
            parent = tree[taxon]
        except KeyError:
            return
        this = taxon

        while True:

            # if shared lineage found
            try:
                idx = lineage.index(this)

                # trim shared lineage
                if idx + 1 < len(lineage):
                    lineage = lineage[slice(idx + 1)]
                break

            # if not found
            except ValueError:
                if parent == this:
                    break
                this = parent
                parent = tree[this]
                continue

    # LCA is the last of shared lineage
    return lineage[-1]
