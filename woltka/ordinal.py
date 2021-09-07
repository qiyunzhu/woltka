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

import numpy as np
from numba import jit
from numba.typed import List

from .align import infer_align_format, assign_parser


def ordinal_mapper(fh, coords, idmap, fmt=None, n=1000000, th=0.8,
                   prefix=False):
    """Read an alignment file and match reads and genes in an ordinal system.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    coords : dict of list
        Gene coordinates table.
    idmap : dict of list
        Gene identifiers.
    fmt : str, optional
        Alignment file format.
    n : int, optional
        Number of lines per chunk.
    th : float
        Minimum threshold of overlap length : alignment length for a match.
    prefix : bool
        Prefix gene IDs with nucleotide IDs.

    See Also
    --------
    align.plain_mapper

    Yields
    ------
    tuple of str
        Query queue.
    dict of set of str
        Subject(s) queue.
    """
    # determine file format
    fmt, head = (fmt, []) if fmt else infer_align_format(fh)

    # assign parser for given format
    parser = assign_parser(fmt)

    # cached list of query Ids for reverse look-up
    # gene Ids are unique, but read Ids can have duplicates (i.e., one read is
    # mapped to multiple loci on a genome), therefore an incremental integer
    # here replaces the original read Id as its identifer
    rids = []
    rid_append = rids.append

    # cached map of read to coordinates
    locmap = defaultdict(list)

    def flush():
        """Match reads in current chunk with genes from all nucleotides.

        Returns
        -------
        tuple of str
            Query queue.
        dict of set of str
            Subject(s) queue.
        """
        # master read-to-gene(s) map
        res = defaultdict(set)

        for nucl, reads in locmap.items():

            # merge and sort coordinates
            # question is to merge an unsorted list into a sorted one
            # it's possible that no gene was annotated on the nucleotide
            try:
                genes = coords[nucl]
            except KeyError:
                continue

            # convert read coordinates to np.array
            reads = np.array(reads, dtype='uint32').reshape(-1, 4)

            # merge gene and read coordinates
            arr = np.concatenate((genes, reads))

            # sort coordinates by 1st column
            # question is to add unsorted read coordinates into already-sorted
            # gene coordinates
            # the 'stable' method uses timsort under the hood, which is
            # efficient for this task
            # TODO: can a np.array be iterated efficiently using sorted
            # indices, without rebuilding a sorted array?
            arr = arr[arr[:, 0].argsort(kind='stable')]

            # get reference to gene identifiers
            gids = idmap[nucl]

            # append prefix if needed
            pfx = nucl + '_' if prefix else ''

            # map reads to genes using the core algorithm
            for read, gene in match_read_gene(arr):

                # merge read-gene pairs to the master map
                res[rids[read]].add(pfx + gids[gene])

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
            locmap = defaultdict(list)

            # next target line number
            target = i + n

        # append read Id, alignment length and location
        idx = len(rids)
        rid_append(query)
        locmap[subject].extend((
            start, -int(-length * th // 1) - 1, 0, idx,
            end, 0, 0, idx))

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
    dict of list of str
        Gene identifiers.
    bool
        Whether there are duplicate gene IDs.

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

    idmap = {}
    gids = None
    gids_append = None

    is_dup = None
    used = set()
    used_add = used.add

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
                gids = idmap[nucl] = []
                gids_append = gids.append
        else:
            x = line.rstrip().split('\t')

            # start and end are based on genome (nucleotide), not gene itself
            try:
                start, end = sorted((int(x[1]), int(x[2])))
            except (IndexError, ValueError):
                raise ValueError(
                    f'Cannot extract coordinates from line: "{line}".')
            gene = x[0]
            idx = len(gids)
            gids_append(gene)
            queue_extend((start, 1, 1, idx,
                          end,   0, 1, idx))

            # check duplicate
            if is_dup is None:
                if idx in used:
                    is_dup = True
                else:
                    used_add(idx)

    # sort gene coordinates per nucleotide
    for nucl, queue in res.items():
        arr = np.array(queue, dtype='uint32').reshape(-1, 4)
        res[nucl] = arr[arr[:, 0].argsort(kind='stable')]
    return res, idmap, is_dup or False


@jit(nopython=True)
def match_read_gene(queue):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : np.array(-1, 4)
        Sorted array of coordinates (loc, is_start, is_gene, idx).

    Yields
    ------
    int
        Read index.
    int
        Gene index.

    See Also
    --------
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

    This function uses Numba's just-in-time compilation to accelerate the
    calculation. Without Numba, its performance will be significantly lower
    (even lower than pure Python without NumPy).

    TODO
    ----
    Cannot cache Numba compilation for unknown reason.
    """
    genes = List()  # current genes
    reads = List()  # current reads

    # cache method references
    genes_append = genes.append
    reads_append = reads.append

    genes_remove = genes.remove
    reads_remove = reads.remove

    # walk through flattened queue of reads and genes
    for loc, is_start, is_gene, idx in queue:
        if is_gene:

            # when a gene starts, add to current genes
            if is_start:
                genes_append((idx, loc))

            # when a gene ends,
            else:

                # find gene start
                for gid, gloc in genes:
                    if gid == idx:
                        genes_remove((gid, gloc))
                        break

                # check current reads
                for rid, rloc, rlen in reads:

                    # is a match if read/gene overlap is long enough
                    if loc - max(gloc, rloc) >= rlen:
                        yield rid, idx

        # the same for reads
        else:
            if is_start:
                reads_append((idx, loc, is_start))
            else:
                for rid, rloc, rlen in reads:
                    if rid == idx:
                        reads_remove((rid, rloc, rlen))
                        break
                for gid, gloc in genes:
                    if loc - max(gloc, rloc) >= rlen:
                        yield idx, gid


def calc_gene_lens(coords, prefix=False):
    """Calculate gene lengths by start and end coordinates.

    Parameters
    ----------
    coords : dict
        Gene coordinates table.
    prefix : bool
        Prefix gene IDs with nucleotide IDs.

    Returns
    -------
    dict of dict
        Mapping of genes to lengths.
    """
    res = {}
    for nucl, queue in coords.items():
        for loc, is_start, _, gid in queue:
            if prefix:
                gid = f'{nucl}_{gid}'
            if is_start:
                res[gid] = 1 - loc
            else:
                res[gid] += loc
    return res
