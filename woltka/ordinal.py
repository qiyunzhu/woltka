#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for matching reads and genes using an ordinal system.
"""

from collections import defaultdict
from itertools import chain
from operator import itemgetter

from .align import infer_align_format, assign_parser


def ordinal_mapper(fh, coords, fmt=None, n=1000000, th=0.8, prefix=False):
    """Read an alignment file and match reads and genes in an ordinal system.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    coords : dict
        Gene coordinates table.
    fmt : str, optional
        Alignment file format.
    n : int, optional
        Number of lines per chunk.
    th : float
        Minimum threshold of overlap length : alignment length for a match.
    prefix : bool
        Prefix gene IDs with Nucleotide IDs.

    See Also
    --------
    align.plain_mapper

    Yields
    ------
    tuple of str
        Query queue.
    dictview of set of str
        Subject(s) queue.
    """
    # determine file format
    fmt, head = (fmt, []) if fmt else infer_align_format(fh)

    # assign parser for given format
    parser = assign_parser(fmt)

    # choose match function depending on whether prefix
    match_func = match_read_gene_pfx if prefix else match_read_gene

    # cached list of query Ids for reverse look-up
    # gene Ids are unique, but read Ids can have duplicates (i.e., one read is
    # mapped to multiple loci on a genome), therefore an incremental integer
    # here replaces the original read Id as its identifer
    rids = []
    rid_append = rids.append

    # cached map of read to alignment length
    lenmap = defaultdict(dict)

    # cached map of read to coordinates
    locmap = defaultdict(list)

    def flush():
        """Match reads in current chunk with genes from all nucleotides.

        Returns
        -------
        tuple of str
            Query queue.
        dictview of set of str
            Subject(s) queue.
        """
        # master read-to-gene(s) map
        res = defaultdict(set)

        for nucl, loci in locmap.items():

            # merge and sort coordinates
            # question is to merge an unsorted list into a sorted one
            # Python's built-in timsort algorithm is efficient at this
            try:
                queue = sorted(coords[nucl] + loci)

            # it's possible that no gene was annotated on the nucleotide
            except KeyError:
                continue

            # map reads to genes using the core algorithm
            for read, gene in match_func(queue, lenmap[nucl], th, nucl):

                # merge read-gene pairs to the master map
                res[rids[read]].add(gene)

        # return matching read Ids and gene Ids
        return res.keys(), res.values()

    this = None  # current query Id
    target = n   # target line number at end of current chunk

    # parse alignment file
    for i, line in enumerate(chain(iter(head), fh)):

        # parse current alignment line
        try:
            query, subject, _, length, start, end = parser(line)[:6]
        except (TypeError, IndexError):
            continue

        # skip if length is not available
        if not length:
            continue

        # when query Id changes and chunk limits has been reached
        if query != this and i >= target:

            # flush: match currently cached reads with genes and yield
            yield flush()

            # re-initiate read Ids, length map and location map
            rids = []
            rid_append = rids.append
            lenmap, locmap = defaultdict(dict), defaultdict(list)

            # next target line number
            target = i + n

        # append read Id, alignment length and location
        idx = len(rids)
        rid_append(query)
        lenmap[subject][idx] = length
        locmap[subject].extend((
            (start, True, False, idx),
            (end,  False, False, idx)))

        this = query

    # final flush
    yield flush()


def ordinal_parser(fh, parser):
    """Alignment parsing functionalities stripped from for `ordinal_mapper`.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    parser : callable
        Function to parse alignment lines of certain format.

    Returns
    -------
    list of str
        Read Ids in same order as in alignment file.
    defaultdict of dict of int
        Map of read indices to alignment lengths per nucleotide.
    defaultdict of list of (int, bool, bool, str)
        Flattened list of read coordinates per nucleotide.

    See Also
    --------
    ordinal_mapper
    read_gene_coords
    .tests.test_ordinal.OrdinalTests.test_ordinal_parser

    Notes
    -----
    This is a dummy function only for test and demonstration purpose but not
    called anywhere in the program. See its unit test for details.

    The data structure is central to the algorithm. See `read_gene_coords` for
    details.
    """
    rids = []
    lenmap = defaultdict(dict)
    locmap = defaultdict(list)

    for line in fh:
        try:
            query, subject, _, length, start, end = parser(line)[:6]
        except (TypeError, IndexError):
            continue
        idx = len(rids)
        rids.append(query)
        lenmap[subject][idx] = length
        locmap[subject].extend((
            (start, True, False, idx),
            (end,  False, False, idx)))

    return rids, lenmap, locmap


def read_gene_coords(fh, sort=False):
    """Read coordinates of genes on genomes.

    Parameters
    ----------
    fh : file handle
        Gene coordinates file.
    sort : bool, optional
        Whether sort gene coordinates.

    Returns
    -------
    dict of list of (int, bool, bool, str)
        Flattened list of gene coordinates per nucleotide.
            Coordinate (nt).
            Whether start (True) or end (False).
            Whether gene (True) or read (False).
            Identifier of gene.

    See Also
    --------
    match_read_gene

    Notes
    -----
    This data structure is central to this algorithm. Starting and ending
    coordinates of each gene are separated and flattened into a sorted list.
    which enables only one round of list traversal for the entire set of genes
    plus reads.
    """
    res = {}
    queue_extend = None
    for line in fh:

        # ">" or "#" indicates genome (nucleotide) name
        c0 = line[0]
        if c0 in '>#':

            # double ">" or "#" indicates genome name, which serves as
            # a super group of subsequent nucleotide names; to be ignored
            if line[1] != c0:
                nucl = line[1:].strip()
                res[nucl] = []
                queue_extend = res[nucl].extend
        else:
            x = line.rstrip().split('\t')

            # start and end are based on genome (nucleotide), not gene itself
            try:
                start, end = sorted((int(x[1]), int(x[2])))
            except (IndexError, ValueError):
                raise ValueError(
                    f'Cannot extract coordinates from line: "{line}".')
            idx = x[0]
            queue_extend(((start, True, True, idx),
                          (end,  False, True, idx)))

    # sort gene coordinates per nucleotide
    if sort:
        sortkey = itemgetter(0)
        for nucl, queue in res.items():
            res[nucl] = sorted(queue, key=sortkey)
    return res


def whether_prefix(coords):
    """Determine whether gene IDs should be prefixed with nucleotide IDs.

    Parameters
    ----------
    coords : dict
        Gene coordinates table.

    Returns
    -------
    bool
        Whether gene IDs should be prefixed.

    See Also
    --------
    read_gene_coords

    Notes
    -----
    It is based on a simple mechanism which checks whether there are duplicate
    gene IDs, and if so, all gene IDs should be prefixed to avoid confusion.
    """
    genes = {}
    for nucl, queue in coords.items():
        for coord, is_start, is_gene, gid in queue:
            if gid not in genes:
                genes[gid] = is_start
            elif genes[gid] == is_start:
                return True
    return False


def match_read_gene(queue, lens, th, pfx=None):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : list of tuple
        Sorted list of elements.
        (loc, is_start, is_gene, id)
    lens : dict
        Read-to-alignment length map.
    th : float
        Threshold for read/gene overlapping fraction.
    pfx : str, optional
        Placeholder for compatibility with match_read_gene_pfx

    Yields
    ------
    int
        Read index.
    str
        Gene ID.

    See Also
    --------
    match_read_gene_pfx
    read_gene_coords
    .tests.test_ordinal.OrdinalTests.test_match_read_gene

    Notes
    -----
    This algorithm is the core of this module. It uses a flattened, sorted
    list to store starting and ending coordinates of both genes and reads.
    Only one round of traversal (O(n)) of this list is needed to accurately
    find all gene-read matches.

    Refer to its unit test `test_match_read_gene` for an actual example and
    illustration.
    """
    genes = {}  # current genes
    reads = {}  # current reads

    # cache method references
    reads_items = reads.items
    genes_items = genes.items

    # walk through flattened queue of reads and genes
    for loc, is_start, is_gene, idx in queue:
        if is_gene:

            # when a gene starts, added to current genes
            if is_start:
                genes[idx] = loc

            # when a gene ends,
            else:

                # check current reads
                for rid, rloc in reads_items():

                    # is a match if read/gene overlap is long enough
                    if loc - max(genes[idx], rloc) + 1 >= lens[rid] * th:
                        yield rid, idx

                # remove it from current genes
                del(genes[idx])

        # the same for reads
        else:
            if is_start:
                reads[idx] = loc
            else:
                for gid, gloc in genes_items():
                    if loc - max(reads[idx], gloc) + 1 >= lens[idx] * th:
                        yield idx, gid
                del(reads[idx])


def match_read_gene_pfx(queue, lens, th, pfx):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : list of tuple
        Sorted list of elements.
        (loc, is_start, is_gene, id)
    lens : dict
        Read-to-alignment length map.
    th : float
        Threshold for read/gene overlapping fraction.
    pfx : str
        Prefix to append to gene IDs.

    Yields
    ------
    int
        Read index.
    str
        Gene ID.

    See Also
    --------
    match_read_gene

    Notes
    -----
    This function is identical to `match_read_gene`, except for that it adds a
    prefix (usually a nucleotide ID) to each gene ID.
    """
    genes, reads = {}, {}
    reads_items, genes_items = reads.items, genes.items
    for loc, is_start, is_gene, idx in queue:
        if is_gene:
            if is_start:
                genes[idx] = loc
            else:
                for rid, rloc in reads_items():
                    if loc - max(genes[idx], rloc) + 1 >= lens[rid] * th:
                        yield rid, f'{pfx}_{idx}'
                del(genes[idx])
        else:
            if is_start:
                reads[idx] = loc
            else:
                for gid, gloc in genes_items():
                    if loc - max(reads[idx], gloc) + 1 >= lens[idx] * th:
                        yield idx, f'{pfx}_{gid}'
                del(reads[idx])
