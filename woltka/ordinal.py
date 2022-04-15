#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Functions for matching reads and genes using an ordinal system.

This modules matches reads and genes based on their coordinates in the genome.
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
    (Paris japonica, 149 billion bp). Therefore it should be sufficient for the
    desired application.

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
from itertools import chain

import numpy as np

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

    # read Ids and effective lengths
    # they are arrays of fixed length = 2 ^ 24
    # there is no need to reset after each flush
    # unlike gene Ids, read Ids can have duplicates (i.e., one read is mapped
    #   to multiple loci in a genome)
    rids = [None] * 2 ** 24
    rels = np.empty(2 ** 24, dtype=np.uint16)

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

        # iterate over nucleotides
        for nucl, locs in locmap.items():

            # it's possible that no gene was annotated on the nucleotide
            try:
                glocs = coords[nucl]
            except KeyError:
                continue

            # get reference to gene identifiers
            gids = idmap[nucl]

            # append prefix if needed
            pfx = nucl + '_' if prefix else ''

            # convert list to array
            locs = np.array(locs, dtype=np.int64)

            # execute ordinal algorithm when reads are many
            # 10 (>5 reads) is an empirically determined cutoff
            if locs.size > 10:

                # map reads to genes using the core algorithm
                for read, gene in match_read_gene(glocs, locs, rels):

                    # add read-gene pairs to the master map
                    res[rids[read]].add(pfx + gids[gene])

            # execute naive algorithm when reads are few
            else:
                for read, gene in match_read_gene_quart(glocs, locs, rels):
                    res[rids[read]].add(pfx + gids[gene])

        # return matching read Ids and gene Ids
        return res.keys(), res.values()

    idx = 0      # current read index
    this = None  # current query Id
    target = n   # target line number at end of current chunk

    # parse alignment file
    for i, line in enumerate(chain(iter(head), fh)):

        # parse current alignment line
        try:
            query, subject, _, length, beg, end = parser(line)[:6]
        except (TypeError, IndexError):
            continue

        # skip if length is not available or zero
        if not length:
            continue

        # when query Id changes and chunk limit has been reached
        if query != this and i >= target:

            # flush: match currently cached reads with genes and yield
            yield flush()

            # reset read index
            idx = 0

            # re-initiate location map
            locmap = defaultdict(list)

            # next target line number
            target = i + n

        # store read Id
        rids[idx] = query

        # store effective length = length * th
        # -int(-x // 1) is equivalent to math.ceil(x) but faster
        rels[idx] = -int(-length * th // 1)

        # add read start and end to queue
        locmap[subject].extend((
            (beg << 24) + idx,
            (end << 24) + (1 << 23) + idx))

        idx += 1
        this = query

    # final flush
    yield flush()


def ordinal_parser_dummy(fh, parser):
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
    match_read_gene_dummy
    .tests.test_ordinal.OrdinalTests.test_ordinal_parser_dummy

    Notes
    -----
    This is a dummy function only for test and demonstration purpose but not
    called anywhere in the program. See its unit test for details.

    The data structure of the queue is:
    0. Coordinate (nt).
    1. Whether start (True) or end (False).
    2. Whether gene (True) or read (False).
    3. Identifier of gene/read.
    """
    rids = []
    rids_append = rids.append
    lenmap = defaultdict(dict)
    locmap = defaultdict(list)

    for line in fh:
        try:
            query, subject, _, length, start, end = parser(line)[:6]
        except (TypeError, IndexError):
            continue
        idx = len(rids)
        rids_append(query)
        lenmap[subject][idx] = length
        locmap[subject].extend((
            (start, True, False, idx),
            (end,  False, False, idx)))

    return rids, lenmap, locmap


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
        Binarized gene coordinate information per genome sequence.
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
                    coords[nucl] = np.array(coords[nucl], dtype=np.int64)
                except KeyError:
                    pass

                # re-initiate gene queue
                nucl = line[1:].strip()
                coords[nucl] = []
                queue_extend = coords[nucl].extend
                gids = idmap[nucl] = []
                gids_append = gids.append

        else:
            x = line.rstrip().split('\t')

            # begin and end positions are based on genome (nucleotide)
            try:
                beg, end = sorted((int(x[1]), int(x[2])))
            except (IndexError, ValueError):
                raise ValueError(
                    f'Cannot extract coordinates from line: "{line}".')
            idx = len(gids)
            gene = x[0]
            gids_append(gene)

            # add positions to queue
            queue_extend(((beg << 24) + (1 << 22) + idx,
                          (end << 24) + (3 << 22) + idx))

            # check duplicate
            if isdup is None:
                if gene in used:
                    isdup = True
                else:
                    used_add(gene)

    # final conversion
    try:
        coords[nucl] = np.array(coords[nucl], dtype=np.int64)
    except KeyError:
        # raise error no data
        pass

    # sort gene coordinates per nucleotide
    if sort:
        for queue in coords.values():
            queue.sort(kind='stable')  # timsort

    return coords, idmap, isdup or False


def match_read_gene(gque, rque, rels):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    gque : np.array(-1, dtype=int64)
        Sorted queue of genes.
    rque : np.array(-1, dtype=int64)
        Paired queue of reads.
    rels : np.array(2 ** 24, dtype=uint16)
        Effective lengths of reads.

    Yields
    ------
    int
        Read index.
    int
        Gene index.

    See Also
    --------
    load_gene_coords
    match_read_gene_dummy
    .tests.test_ordinal.OrdinalTests.test_match_read_gene_dummy

    Notes
    -----
    This algorithm is the core of this module. It uses a flattened, sorted
    list to store starting and ending coordinates of both genes and reads.
    Only one round of iteration (O(n)) over this list is needed to accurately
    find all gene-read matches.

    Refer to its unit test `test_match_read_gene` for an illustrated example.

    This function is the most compute-intensive step in the entire analysis,
    therefore it has been deeply optimized to increase performance wherever
    possible.

    Note: Repeated bitwise operations are usually more efficient that a single
    bitwise operation assigned to a new variable.
    """

    # merge pre-sorted genes with reads with unknown sorting status
    # timsort is efficient for this task
    queue = np.sort(np.concatenate([gque, rque]), kind='stable')

    genes = {}  # current genes cache
    reads = {}  # current reads cache

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
                gloc = genes_pop(code & (1 << 22) - 1)

                # check cached reads for matches
                for rid, rloc in reads_items():

                    # is a match if read/gene overlap is long enough
                    if (code >> 24) - max(gloc, rloc) >= rels[rid] - 1:
                        yield rid, code & (1 << 22) - 1

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
                # `dropwhile` unmatched, then yield all remaining ones
                for gid, gloc in genes_items():

                    # same as above
                    if (code >> 24) - max(gloc, rloc) >= rels[rid] - 1:
                        yield rid, gid


def match_read_gene_naive(gque, rque, rels):
    """Associate reads with genes using a naive approach, which performs
    nested iteration over genes and reads.

    Parameters
    ----------
    gque : np.array(-1, dtype=int64)
        Sorted queue of genes.
    rque : np.array(-1, dtype=int64)
        Paired queue of reads.
    rels : np.array(2 ** 24, dtype=uint16)
        Effective lengths of reads.

    Yields
    ------
    int
        Read index.
    int
        Gene index.

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
    # only genes are cached
    genes = {}
    genes_pop = genes.pop

    # iterate over reads (paired)
    it = iter(rque)
    for x, y in zip(it, it):

        # pre-calculate read metrics
        rid = x & (1 << 22) - 1  # index
        beg = x >> 24            # start coordinate
        end = y >> 24            # end coordinate
        L = rels[rid] - 1        # effective length - 1

        # iterate over gene positions
        for code in gque:

            # at start, add gene to cache
            if not code & (1 << 23):
                genes[code & (1 << 22) - 1] = code >> 24

            # at end, check overlap while removing gene from cache:
            #   min(gene end, read end) - max(gene start, read start) >=
            #   effective length - 1
            elif (min(code >> 24, end) -
                  max(genes_pop(code & (1 << 22) - 1), beg)) >= L:
                yield rid, code & (1 << 22) - 1


def match_read_gene_quart(gque, rque, rels):
    """Associate reads with genes by iterating between read region and nearer
    end of gene queue.

    Parameters
    ----------
    gque : np.array(-1, dtype=int64)
        Sorted queue of genes.
    rque : np.array(-1, dtype=int64)
        Paired queue of reads.
    rels : np.array(2 ** 24, dtype=uint16)
        Effective lengths of reads.

    Yields
    ------
    int
        Read index.
    int
        Gene index.

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
    over half of the genome, and they cannot be detected within out iterating
    to end of the gene queue.

    This method uses bisection to find insertion points of read start in the
    pre-sorted gene queue, which has O(logn) time.
    """
    n = gque.size  # entire search space
    mid = n // 2   # mid point

    # iterate over paired read starts and ends
    it = iter(rque)
    for x, y in zip(it, it):

        # pre-calculate read metrics
        rid = x & (1 << 22) - 1  # index
        beg = x >> 24            # start coordinate
        end = y >> 24            # end coordinate
        L = rels[rid] - 1        # effective length - 1

        # genes starting within read region
        # will record their coordinates
        within = {}
        within_pop = within.pop

        # locate read start using bisection
        i = gque.searchsorted(x, side='right')

        # read starts in left half of gene queue
        if i < mid:

            # genes starting before read region
            # no need to record coordinates as overlap is not possible
            before = set()
            before_add = before.add
            before_remove = before.remove

            # iterate over gene positions before read region
            for code in gque[:i]:

                # add gene to cache when it starts
                if not code & (1 << 23):
                    before_add(code & (1 << 22) - 1)

                # drop gene from cache when it ends
                else:
                    before_remove(code & (1 << 22) - 1)

            # iterate over gene positions after read start
            for code in gque[i:]:

                # stop when exceeding read end
                if code > y:
                    break

                # when gene starts, add it and its coordinate to cache
                elif not code & (1 << 23):
                    within[code & (1 << 22) - 1] = code >> 24

                # when gene ends, check overlap and remove it from cache
                else:
                    gid = code & (1 << 22) - 1

                    # most genes should start before read region
                    try:
                        before_remove(gid)

                    # if gene started within read region,
                    # overlap = gene end - gene start
                    except KeyError:
                        if (code >> 24) - within_pop(gid) >= L:
                            yield rid, gid

                    # if gene started before read region,
                    # overlap = gene end - read start
                    else:
                        if (code >> 24) - beg >= L:
                            yield rid, gid

            # yield all genes that started before read but not yet ended
            for gid in before:
                yield rid, gid

            # yield genes that started after read but not yet ended,
            # until overlap = read end - gene start is too short
            for gid, gloc in within.items():
                if end - gloc < L:
                    break
                yield rid, gid

        # read starts in right half of gene queue
        else:

            # genes starting after read region
            # no need to record coordinates as overlap is not possible
            after = set()
            after_add = after.add
            after_remove = after.remove

            # iterate over gene positions within read region
            while i < n:
                code = gque[i]

                # stop when exceeding read end
                if code > y:
                    break

                # when gene starts, add it and its coordinate to cache
                if not code & (1 << 23):
                    within[code & (1 << 22) - 1] = code >> 24

                # when gene ends, check overlap and remove it from cache
                # gene end must be <= read end
                # if gene not found in cache (meaning gene started before read
                # region), use read start, otherwise use gene start
                elif (code >> 24) - within_pop(code & (1 << 22) - 1, beg) >= L:
                    yield rid, code & (1 << 22) - 1

                # move to next position
                i += 1

            # iterate over gene positions after read region
            for code in gque[i:]:

                # add gene to cache when it starts
                if not code & (1 << 23):
                    after_add(code & (1 << 22) - 1)

                # when gene ends,
                else:

                    # most remaining genes should start after read region
                    try:
                        after_remove(code & (1 << 22) - 1)

                    # otherwise, the gene spans over entire read region
                    except KeyError:

                        # check overlap while removing gene from cache
                        # gene end must be >= read end
                        # if gene not found in cache, use read start, otherwise
                        # use gene start
                        if end - within_pop(code & (1 << 22) - 1, beg) >= L:
                            yield rid, code & (1 << 22) - 1


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
    for loc, is_start, is_gene, idx in queue:
        if is_gene:

            # when a gene starts, added to gene cache
            if is_start:
                genes[idx] = loc

            # when a gene ends,
            else:

                # find gene start and remove it from cache
                gloc = genes.pop(idx)

                # check cached reads for matches
                for rid, rloc in reads.items():

                    # is a match if read/gene overlap is long enough
                    if loc - max(gloc, rloc) + 1 >= lens[rid] * th:
                        yield rid, idx

        # the same for reads
        else:
            if is_start:
                reads[idx] = loc
            else:
                rloc = reads.pop(idx)
                for gid, gloc in genes.items():
                    if loc - max(rloc, gloc) + 1 >= lens[idx] * th:
                        yield idx, gid


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
            # length = end - start + 1
            if not code & (1 << 23):  # start
                res[gid] = 1 - (code >> 24)
            else:  # end
                res[gid] += code >> 24
    return res
