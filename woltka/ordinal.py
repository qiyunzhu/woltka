#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for matching reads and genes using an ordinal system.

This module matches reads and genes based on their coordinates in the genome.
It uses an efficient algorithm, in which starting and ending positions of reads
and genes are flattened and sorted by coordinate into one sequence, followed by
one iteration over the sequence to identify all read-gene overlaps.

This algorithm was inspired by the sweep line algorithm [1] in computational
geometry.

    [1] Shamos, Michael Ian, and Dan Hoey. "Geometric intersection problems."
        17th Annual Symposium on Foundations of Computer Science (sfcs 1976).
        IEEE, 1976.

In this sequence, each position contains 4 pieces of information:

    0. Coordinate (nt).
    1. Whether start or end.
    2. Whether gene or read.
    3. Index of gene / read.

For efficient computing, the information is stored as a 64-bit signed integer
(int64), and accessed using bitwise operations. The bits are defined as (from
right to left):

    Bits  1-22 (22): Index of gene / read (max.: 4,194,303).
    Bits    23  (1): Whether it is a read (0) or a gene (1).
    Bits    24  (1): Whether it is a start (0) or an end (1).
    Bits 25-63 (39): Coordinate (max.: 549,755,813,887).

To construct a position:

    code = (nt << 24) + is_gene * (1 << 23) + is_end * (1 << 22) + idx

To extract information from a position:

    Index:        code & (1 << 22) - 1
    Gene or read: code & (1 << 22)
    End or start: code & (1 << 23)
    Coordinate:   code >> 24

Rationales of this design:

    Index (max. = 4 million): It is larger than the largest number of protein-
    coding genes in a genome (Trichomonas vaginalis, ca. 60,000). On the other
    side, the input chunk size should be no more than 4 million, otherwise it
    is possible that reads mapped to a single genome exceed the upper limit of
    indices.

    Read (0) or gene (1): Adding one bit costs some compute. It should be done
    during gene queue construction, which only happens during database loading,
    instead of during read queue construction, which happens for every chunk of
    input data.

    Start (0) or end (1): This ensures that the points (uint64) can be directly
    sorted in ascending order, such that start always occurs before end, even
    if their coordinates are equal (i.e., the gene/read has a length of 1 bp).

    Coordinate (max. = 550 billion): It is larger than the largest known genome
    (Tmesipteris oblanceolata, 160 Gbp). Therefore it should be sufficient for
    the desired application.

Notes:

    If data type is uint64 instead of int64, the maximum coordinate can be 2x
    as large. However, NumPy does not allow bitwise operations on uint64.

    To remove the upper limit of coordinate, one may remove `dtype=np.int64`
    from the code. This will slightly slow down the code even if no coordinate
    exceeds the upper limit, and it will notably reduce efficiency when some
    coordinates exceed the upper limit, because the array will downgrade to
    data type to 'O' (Python object), which is inefficient.
"""

from collections import defaultdict

import numpy as np
from numba import njit
from numba.typed import Dict
from numba.types import uint16, int64, boolean
from numba import config

from .align import iter_align

nojit = config.DISABLE_JIT


def match_read_gene_dummy(queue, lens, th):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : list of tuple of (int, bool, bool, int)
        Sorted queue of coordinates (location, start or end, gene or read,
        index).
    lens : dict
        Read-to-alignment length map.
    th : float
        Threshold for read/gene overlapping fraction.

    Yields
    ------
    int
        Read index.
    int
        Gene index.

    See Also
    --------
    match_read_gene
    .tests.test_ordinal.OrdinalTests.test_match_read_gene_dummy

    Notes
    -----
    This is a dummy function which is only for test and demonstration purpose,
    but is not called anywhere in the program.

    The formal function `match_read_gene` extensively uses bitwise operations
    and thus is hard to read. Therefore the current function, which represents
    the original form of prior to optimization, is retained.
    """
    genes = {}  # current genes
    reads = {}  # current reads

    # walk through flattened queue of reads and genes
    for loc, is_end, is_gene, idx in queue:
        if is_gene:

            # when a gene starts, added to gene cache
            if not is_end:
                genes[idx] = loc

            # when a gene ends,
            else:

                # find gene start and remove it from cache
                gloc = genes.pop(idx)

                # check cached reads for matches
                for rid, rloc in reads.items():

                    # is a match if read/gene overlap is long enough
                    if loc - max(gloc, rloc) >= lens[rid] * th:
                        yield rid, idx

        # the same for reads
        else:
            if not is_end:
                reads[idx] = loc
            else:
                rloc = reads.pop(idx)
                for gid, gloc in genes.items():
                    if loc - max(rloc, gloc) >= lens[idx] * th:
                        yield idx, gid


def ordinal_mapper(fh, coords, idmap, fmt=None, excl=None, n=2**20, th=0.8,
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
    excl : set, optional
        Subjects to exclude.
    n : int, optional
        Number of reads per chunk (max. 2 ** 22 = 4,194,304).
    th : float
        Minimum threshold of overlap length : alignment length for a match.
    prefix : bool
        Prefix gene IDs with nucleotide IDs.

    See Also
    --------
    .align.plain_mapper

    Yields
    ------
    list of str
        Query queue.
    list of set of str
        Subject(s) queue.
    """
    it = iter_align(fh, fmt, excl, True)

    # arguments for flush_chunk
    args = (coords, idmap, th, prefix)

    # cached lists of read Ids and lengths (pre-allocate space)
    # gene Ids are unique, but read Ids can have duplicates (i.e., one read is
    # mapped to multiple loci on a genome), therefore an incremental integer
    # here replaces the original read Id as its identifer
    rids = [None] * n
    lens = np.empty((n,), dtype=np.uint16)

    # cached map of reads to per-genome coordinates
    locmap = defaultdict(list)

    # current read index in the cached lists; will reset after each flush
    idx = 0

    # parse alignment file
    for query, records in it:

        # exclude hits with unavailable or zero length
        records = [x for x in records if x[2]]

        # when chunk limit is about to be exceeded by the next query, match
        # currently cached reads with genes, flush, and reset
        if idx + len(records) > n:
            yield flush_chunk(idx, locmap, rids, lens, *args)
            locmap = defaultdict(list)
            idx = 0

        # extract alignment info and add to cache
        for subject, _, length, beg, end in records:
            rids[idx] = query
            lens[idx] = length
            locmap[subject].extend((
                (beg << 24) + idx,
                (end << 24) + (1 << 23) + idx))
            idx += 1

    # final flush
    yield flush_chunk(idx, locmap, rids, lens, *args)


def flush_chunk(n, rlocmap, rids, rlens, glocmap, gidmap, th, prefix):
    """Match reads in current chunk with genes from all genomes.

    Parameters
    ----------
    n : int
        Number of reads to flush.
    rlocmap : dict of list
        Read coordinates per genome.
    rids : list of str
        Read IDs.
    rlens : 1-D array_like of int
        Read lengths.
    glocmap : dict of list
        Gene coordinates per genome.
    gidmap : dict of list
        Gene identifiers.
    th : float
        Length threshold.
    prefix : bool
        Prefix gene IDs with nucleotide IDs.

    Returns
    -------
    list of str
        Query queue.
    list of set of str
        Subject(s) queue.
    """
    # master read-to-gene(s) map
    res = defaultdict(set)

    # effective length = length * th
    rels = np.ceil(rlens[:n] * th).astype(np.uint16)

    # iterate over nucleotides
    for nucl, rlocs in rlocmap.items():

        # it's possible that no gene was annotated on the nucleotide
        try:
            glocs = glocmap[nucl]
        except KeyError:
            continue

        # get reference to gene identifiers
        gids = gidmap[nucl]

        # append prefix if needed
        pfx = nucl + '_' if prefix else ''

        # convert list to array
        rlocs = np.array(rlocs, dtype=np.int64)

        # execute ordinal algorithm when reads are many
        # 10 (>5 reads) is an empirically determined cutoff
        if rlocs.size > 10:

            # merge pre-sorted genes with reads of unknown sorting status
            queue = np.concatenate((glocs, rlocs))

            # sort genes and reads into a mixture
            # timsort is efficient for this task
            queue.sort(kind='stable')

            # map reads to genes using the core algorithm
            matches = match_read_gene(queue, rels)

        # execute naive algorithm when reads are few
        else:
            matches = match_read_gene_quart(glocs, rlocs, rels)

        # add read-gene pairs to the master map
        for read, gene in matches:
            res[rids[read]].add(pfx + gids[gene])

    # return matching read Ids and gene Ids
    return res.keys(), res.values()


def load_gene_coords(fh, sort=False):
    """Read coordinates of genes on genomes.

    Parameters
    ----------
    fh : file handle
        Gene coordinates file.
    sort : bool, optional
        Whether sort gene coordinates.

    Returns
    -------
    dict of np.array(-1, dtype=int64)
        Binarized gene position information per genome.
    dict of list of str
        Gene IDs.
    bool
        Whether there are duplicate gene IDs.

    See Also
    --------
    match_read_gene
    """
    coords = {}
    queue_extend = None

    idmap = {}
    gids = None
    gids_append = None

    isdup = None
    used = set()
    used_add = used.add

    nucl = None  # current nucleotide

    for line in fh:

        # ">" or "#" indicates genome (nucleotide) name
        c0 = line[0]
        if c0 in '>#':

            # double ">" or "#" indicates genome name, which serves as
            # a super group of subsequent nucleotide names; to be ignored
            if line[1] != c0:

                # convert gene queue to np.array
                try:
                    coords[nucl] = encode_genes(coords[nucl])
                except KeyError:
                    pass

                # re-initiate gene queue
                nucl = line[1:].strip()
                coords[nucl] = []
                queue_extend = coords[nucl].extend
                gids = idmap[nucl] = []
                gids_append = gids.append

        else:

            # extract Id, start, end
            try:
                gene, beg, end = line.rstrip().split('\t')
            except ValueError:
                raise ValueError(
                    f'Cannot extract coordinates from line: "{line}".')

            # add positions to queue
            queue_extend((beg, end))

            # record gene Id
            gids_append(gene)

            # check duplicate
            if isdup is None:
                if gene in used:
                    isdup = True
                else:
                    used_add(gene)

    # final conversion
    try:
        coords[nucl] = encode_genes(coords[nucl])
    except KeyError:
        raise ValueError('No coordinate was read from file.')

    # sort gene coordinates per nucleotide
    if sort:
        for queue in coords.values():
            queue.sort(kind='stable')  # timsort

    return coords, idmap, isdup or False


def encode_genes(lst):
    """Encode gene positions into a binary queue.

    Parameters
    ----------
    lst : list of int
        Flattened list of start and end coordinates.

    Returns
    -------
    np.array(-1, dtype=int64)
        Encoded gene queue.
    """
    try:
        arr = np.asarray(lst, dtype=np.int64)
    except ValueError:
        raise ValueError('Invalid coordinate(s) found.')

    n = arr.size // 2
    idx = np.arange(n)  # gene indices

    # separate start (odd) and end (even) positions
    beg, end = arr[0::2], arr[1::2]

    # order each pair of start and end coordinates such that smaller one
    # comes first
    # faster than np.sort since there are only two numbers
    # < is slightly faster than np.less
    cmp = beg < end
    lo = np.where(cmp, beg, end)
    hi = np.where(cmp, end, beg)

    # encode coordinate, start/end, is gene, and index into one integer
    lo = np.left_shift(lo - 1, 24) + (1 << 22) + idx
    hi = np.left_shift(hi, 24) + (3 << 22) + idx

    # fastest way to interleave two arrays
    # https://stackoverflow.com/questions/5347065/
    que = np.empty((2 * n,), dtype=np.int64)
    que[0::2] = lo
    que[1::2] = hi

    return que


@njit((int64[:], int64[:], uint16[:]))
def match_read_gene(queue, rels):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : np.array(-1, dtype=int64)
        Paired queue of reads.
    rels : np.array(-1, dtype=uint16)
        Effective lengths of reads.

    Returns
    ------
    list of tuple of (int, int)
        Indices of matched reads and genes.

    See Also
    --------
    load_gene_coords
    match_read_gene_dummy
    .tests.test_ordinal.OrdinalTests.test_match_read_gene_dummy

    Notes
    -----
    This algorithm is the core of this module. It uses a flattened, sorted
    list to store starting and ending coordinates of both genes and reads.
    Only one round of traversal (O(n)) of this list is needed to accurately
    find all gene-read matches.

    Refer to its unit test `test_match_read_gene` for an illustrated example.

    This function is the most compute-intensive step in the entire analysis,
    therefore it has been deeply optimized to increase performance wherever
    possible.

    Note: Repeated bitwise operations are usually more efficient that a single
    bitwise operation assigned to a new variable.
    """
    # matched read-gene pairs
    # numba has issues with jitting generator functions
    # see https://github.com/numba/numba/issues/6993
    # therefore this function returns rather than yields
    matches = []
    matches_append = matches.append

    # caches of current genes and reads
    genes = {} if nojit else Dict.empty(key_type=int64, value_type=int64)
    reads = {} if nojit else Dict.empty(key_type=int64, value_type=int64)

    # cache method references
    genes_items = genes.items
    reads_items = reads.items

    genes_pop = genes.pop
    reads_pop = reads.pop

    # walk through flattened queue of reads and genes
    for code in queue:

        # if this is a gene,
        if code & (1 << 22):

            # when a gene starts,
            if not code & (1 << 23):

                # add it to cache
                genes[code & (1 << 22) - 1] = code >> 24

            # when a gene ends,
            else:

                # find the gene's start
                gid = code & (1 << 22) - 1
                gloc = genes_pop(gid)

                # check cached reads for matches
                for rid, rloc in reads_items():

                    # is a match if read/gene overlap is long enough
                    if (code >> 24) - max(gloc, rloc) >= rels[rid]:
                        matches_append((rid, gid))

        # if this is a read,
        else:

            # when a read starts,
            if not code & (1 << 23):

                # add it to cache
                reads[code & (1 << 22) - 1] = code >> 24

            # when a read ends,
            else:

                # find the read's start
                rid = code & (1 << 22) - 1
                rloc = reads_pop(rid)

                # check cached genes
                # theoretically, an optimization is to loop in reverse order,
                # `dropwhile` unmatched, then return all remaining ones
                for gid, gloc in genes_items():

                    # same as above
                    if (code >> 24) - max(gloc, rloc) >= rels[rid]:
                        matches_append((rid, gid))
    return matches


@njit((int64[:], int64[:], uint16[:]))
def match_read_gene_naive(gque, rque, rels):
    """Associate reads with genes using a naive approach, which performs
    nested iteration over genes and reads.

    Parameters
    ----------
    gque : np.array(-1, dtype=int64)
        Sorted queue of genes.
    rque : np.array(-1, dtype=int64)
        Paired queue of reads.
    rels : np.array(-1, dtype=uint16)
        Effective lengths of reads.

    Returns
    ------
    list of tuple of (int, int)
        Indices of matched reads and genes.

    See Also
    --------
    match_read_gene

    Notes
    -----
    This is a reference implementation of the naive nested iteration method
    with limited optimization. It is O(nm), where n and m are the numbers of
    genes and reads, respectively. It should be much slower than the ordinal
    method. However, when the number of reads is small, it may be efficient
    because it saves the sorting step.
    """
    matches = []
    matches_append = matches.append

    # only genes are cached
    genes = {} if nojit else Dict.empty(key_type=int64, value_type=int64)
    genes_pop = genes.pop

    # iterate over reads in pairs
    # this way is faster than rque.reshape(-1, 2)
    it = iter(rque)
    for x, y in zip(it, it):

        # pre-calculate read metrics
        rid = x & (1 << 22) - 1  # index
        beg = x >> 24            # start coordinate
        end = y >> 24            # end coordinate
        L = rels[rid]            # effective length

        # iterate over gene positions
        for code in gque:

            # at start, add gene to cache
            if not code & (1 << 23):
                genes[code & (1 << 22) - 1] = code >> 24

            # at end, check overlap while removing gene from cache:
            #   min(gene end, read end) - max(gene start, read start) >=
            #   effective length
            elif (min(code >> 24, end) -
                  max(genes_pop(code & (1 << 22) - 1), beg)) >= L:
                matches_append((rid, code & (1 << 22) - 1))
    return matches


@njit((int64[:], int64[:], uint16[:]))
def match_read_gene_quart(gque, rque, rels):
    """Associate reads with genes by iterating between read region and nearer
    end of gene queue.

    Parameters
    ----------
    gque : np.array(-1, dtype=int64)
        Sorted queue of genes.
    rque : np.array(-1, dtype=int64)
        Paired queue of reads.
    rels : np.array(-1, dtype=uint16)
        Effective lengths of reads.

    Returns
    ------
    list of tuple of (int, int)
        Indices of matched reads and genes.

    See Also
    --------
    match_read_gene_naive

    Notes
    -----
    This optimized method is similar to the naive method, but it reduces the
    search space from n to n / 4, thus it significantly improves efficiency.

    The idea is that one only needs to iterate between the position where read
    can fit it, and the nearer end of the gene queue.

    Further reduction of search space is not possible, because genes may span
    over half of the genome, and they cannot be detected without iterating to
    end of the gene queue.

    This method uses bisection to find the insertion point of read start in the
    pre-sorted gene queue, which has O(logn) time.
    """
    matches = []
    matches_append = matches.append

    # genes within read region
    within = {} if nojit else Dict.empty(key_type=int64, value_type=int64)
    within_pop = within.pop

    # genes outside read region
    # there is some numba issue, otherwise set may be more suitable
    outside = {} if nojit else Dict.empty(key_type=int64, value_type=boolean)
    outside_pop = outside.pop

    n = gque.size  # entire search space
    mid = n // 2   # mid point

    # iterate over paired read starts and ends
    it = iter(rque)
    for x, y in zip(it, it):

        # pre-calculate read metrics
        rid = x & (1 << 22) - 1  # index
        beg = x >> 24            # start coordinate
        end = y >> 24            # end coordinate
        L = rels[rid]            # effective length

        # locate read start using bisection
        i = np.searchsorted(gque, x, side='right')

        # read starts in left half of gene queue
        if i < mid:

            # iterate over gene positions before read region
            for code in gque[:i]:

                # add gene to cache when it starts
                if not code & (1 << 23):
                    outside[code & (1 << 22) - 1] = True

                # drop gene from cache when it ends
                else:
                    outside_pop(code & (1 << 22) - 1)

            # iterate over gene positions after read start
            for code in gque[i:]:

                # stop when exceeding read end
                if code > y:
                    break

                # when gene starts, add it and its coordinate to cache
                elif not code & (1 << 23):
                    within[code & (1 << 22) - 1] = code >> 24

                # when gene ends, check overlap and remove it from cache
                # if gene started before read region (most cases),
                # overlap = gene end - read start
                elif outside_pop(code & (1 << 22) - 1, False):
                    if (code >> 24) - beg >= L:
                        matches_append((rid, code & (1 << 22) - 1))

                # if gene started within read region,
                # overlap = gene end - gene start
                elif (code >> 24) - within_pop(code & (1 << 22) - 1) >= L:
                    matches_append((rid, code & (1 << 22) - 1))

            # return all genes that started before read but not yet ended
            for gid in outside:
                matches_append((rid, gid))
            outside.clear()

            # return genes that started after read but not yet ended,
            # until overlap = read end - gene start is too short
            for gid, gloc in within.items():
                if end - gloc < L:
                    break
                matches_append((rid, gid))
            within.clear()

        # read starts in right half of gene queue
        else:

            # iterate over gene positions within read region
            while i < n:
                code = gque[i]

                # stop when exceeding read end
                if code > y:
                    break
                gid = code & (1 << 22) - 1

                # when gene starts, add it and its coordinate to cache
                if not code & (1 << 23):
                    within[gid] = code >> 24

                # when gene ends, check overlap and remove it from cache
                # gene end must be <= read end
                # if gene not found in cache (gene started before read region),
                # use read start, otherwise use gene start
                elif (code >> 24) - within_pop(gid, beg) >= L:
                    matches_append((rid, gid))

                # move to next position
                i += 1

            # iterate over gene positions after read region
            for code in gque[i:]:
                gid = code & (1 << 22) - 1

                # add gene to cache when it starts
                if not code & (1 << 23):
                    outside[gid] = True

                # when gene ends,
                # most remaining genes should start after read region
                # otherwise,
                elif not outside_pop(gid, False):

                    # check overlap while removing gene from cache
                    # gene end must be >= read end
                    # if gene not found in cache, use read start, otherwise
                    # use gene start
                    if end - within_pop(gid, beg) >= L:
                        matches_append((rid, gid))
    return matches


def calc_gene_lens(mapper):
    """Calculate gene lengths by start and end coordinates.

    Parameters
    ----------
    mapper : callable
        Ordinal mapper.

    Returns
    -------
    dict of int
        Mapping of genes to lengths.
    """
    res = {}
    prefix = mapper.keywords['prefix']
    idmap = mapper.keywords['idmap']
    for nucl, queue in mapper.keywords['coords'].items():
        idmap_ = idmap[nucl]
        nucl += '_'
        for code in queue:
            gid = idmap_[code & (1 << 22) - 1]
            if prefix:
                gid = nucl + gid
            if not code & (1 << 23):  # start
                res[gid] = -(code >> 24)
            else:  # end
                res[gid] += code >> 24
    return res


def load_gene_lens(fh):
    """Directly load gene lengths from a coordinates file.

    Parameters
    ----------
    fh : file handle
        Gene coordinates file.

    Returns
    -------
    dict of int
        Gene ID to length mapping.

    See Also
    --------
    load_gene_coords

    Notes
    -----
    Used by the "normalize" command. This is more efficient than:

    >>> coords, idmap, prefix = load_gene_coords(f, sort=True)
    >>> mapper = partial(ordinal_mapper, coords=coords, idmap=idmap,
                         prefix=prefix)
    >>> sizemap = calc_gene_lens(mapper)

    """
    lens, lens_ = {}, None
    used, isdup = set(), None
    used_add = used.add
    for line in fh:
        c0 = line[0]
        if c0 in '>#':
            if line[1] != c0:
                nucl = line[1:].strip()
                lens_ = lens[nucl] = []
        else:
            x = line.rstrip().split('\t')
            try:
                gene, beg, end = x[0], int(x[1]), int(x[2])
            except (IndexError, ValueError):
                raise ValueError(
                    f'Cannot calculate gene length from line: "{line}".')
            lens_.append((gene, abs(end - beg) + 1))
            if isdup is None:
                if gene in used:
                    isdup = True
                else:
                    used_add(gene)
    if isdup:
        return {f'{k}_{g}': L for k, v in lens.items() for g, L in v}
    else:
        return dict(e for k, v in lens.items() for e in v)
