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
from bisect import bisect

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
    parser = assign_parser(fmt, ext=True)

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

            # execute ordinal algorithm when reads are many
            # 8 (5+ reads) is an empirically determined cutoff
            if len(locs) > 8:

                # merge and sort coordinates
                # question is to add unsorted read coordinates into pre-sorted
                # gene coordinates
                # Python's Timsort algorithm is efficient for this task
                queue = sorted(chain(glocs, locs))

                # map reads to genes using the core algorithm
                for read, gene in match_read_gene(queue):

                    # add read-gene pairs to the master map
                    res[rids[read]].add(pfx + gids[gene])

            # execute naive algorithm when reads are few
            else:
                for read, gene in match_read_gene_quart(glocs, locs):
                    res[rids[read]].add(pfx + gids[gene])

        # return matching read Ids and gene Ids
        return res.keys(), res.values()

    this = None  # current query Id
    target = n   # target line number at end of current chunk

    # parse alignment file
    for i, row in enumerate(parser(chain(iter(head), fh))):
        query, subject, _, length, beg, end = row[:6]

        # skip if length is not available or zero
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

        # effective length = length * th
        # -int(-x // 1) is equivalent to math.ceil(x) but faster
        # this value must be >= 1
        locmap[subject].extend((
            (beg << 48) + (-int(-length * th // 1) << 31) + idx,
            (end << 48) + idx))
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
    """
    rids = []
    lenmap = defaultdict(dict)
    locmap = defaultdict(list)

    for row in parser(fh):
        query, subject, _, length, start, end = row[:6]
        idx = len(rids)
        rids.append(query)
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
    dict of int
        Binarized gene coordinate information per nucleotide.
    dict of list of str
        Gene IDs.
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

    See the docstring of `match_read_gene` for details.
    """
    coords = {}
    queue_extend = None

    idmap = {}
    gids = None
    gids_append = None

    isdup = None
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
                coords[nucl] = []
                queue_extend = coords[nucl].extend
                gids = idmap[nucl] = []
                gids_append = gids.append
        else:
            x = line.rstrip().split('\t')

            # begin and end positions are based on genome (nucleotide)
            try:
                beg, end = int(x[1]), int(x[2])
            except (IndexError, ValueError):
                raise ValueError(
                    f'Cannot extract coordinates from line: "{line}".')
            idx = len(gids)
            gene = x[0]
            gids_append(gene)
            if beg > end:
                beg, end = end, beg
            queue_extend(((beg << 48) + (3 << 30) + idx,
                          (end << 48) + (1 << 30) + idx))

            # check duplicate
            if isdup is None:
                if gene in used:
                    isdup = True
                else:
                    used_add(gene)

    # sort gene coordinates per nucleotide
    if sort:
        for queue in coords.values():
            queue.sort()

    return coords, idmap, isdup or False


def match_read_gene(queue):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : list of int
        Sorted queue of coordinates.

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
    Only one round of traversal (O(n)) of this list is needed to accurately
    find all gene-read matches.

    Refer to its unit test `test_match_read_gene` for an illustrated example.

    This function is the most compute-intensive step in the entire analysis,
    therefore it has been deeply optimized to increase performance wherever
    possible. Notably, it extensively uses bitwise operations to extract
    multiple pieces of information from a single integer.

    Specifically, each coordinate (an integer) has the following information
    (from right to left):

    - Bits  1-30: Index of gene / read (30 bits, max.: 1,073,741,823).
    - Bits    31: Whether it is a gene (1) or a read (0) (1 bit).
    - Bits 32-58: Whether it is the start (positive) or end (0) of a gene /
                  read. If start, the value represents the effective length of
                  an alignment if it's a read, or 1 if it's a gene (17 bits,
                  max.: 131,071).
    - Bits 59-  : Coordinate (position on the genome, nt) (unlimited)

    The Python code to extract these pieces of information is:

    - Coordinate:       `code >> 48`
    - Effective length: `code >> 31 & (1 << 17) - 1`, or `code >> 31 & 131071`
    - Gene or read:     `code & (1 << 30)`
    - Gene/read index:  `code & (1 << 30) - 1`

    Note: Repeated bitwise operations are usually more efficient that a single
    bitwise operation assigned to a new variable.
    """
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
        if code & (1 << 30):

            # when a gene begins,
            # if code >> 31 & 131071:
            if code & (1 << 31):

                # add its index and coordinate to cache
                genes[code & (1 << 30) - 1] = code >> 48

            # when a gene ends,
            else:

                # find gene start
                gloc = genes_pop(code & (1 << 30) - 1)

                # check cached reads for matches
                for rid, rloc in reads_items():

                    # is a match if read/gene overlap is long enough
                    #   code >> 48:     read end coordinate
                    #   gloc:           gene start coordinate
                    #   rloc >> 17`:    read start coordinate
                    #   rloc & 131071`: effective length - 1
                    if (code >> 48) - max(gloc, rloc >> 17) >= rloc & 131071:
                        yield rid, code & (1 << 30) - 1

        # if this is a read,
        else:

            # when a read begins,
            if code >> 31 & 131071:

                # add its index, coordinate and effective length - 1 to cache
                # the latter two are stored as a single integer
                reads[code & (1 << 30) - 1] = (code >> 31) - 1

            # when a read ends,
            else:

                # find read start and effective length
                rloc = reads_pop(code & (1 << 30) - 1)

                # check cached genes
                # a potential optimization is to pre-calculate `rloc >> 17`
                #   and `(code >> 48) - (rloc & 131071)`, however, it is not
                #   worth the overhead in real applications because there is
                #   usually zero or one gene in cache
                # another potential optimization is to replace `max` (which is
                #   a function call with a ternary operator, but one needs the
                #   first optimization prior to this
                for gid, gloc in genes_items():

                    # same as above
                    if (code >> 48) - max(gloc, rloc >> 17) >= rloc & 131071:
                        yield code & (1 << 30) - 1, gid


def match_read_gene_naive(geneque, readque):
    """Associate reads with genes using a naive approach, which performs
    nested iteration over genes and reads.

    Parameters
    ----------
    geneque : list of int
        Sorted queue of genes.
    readque : list of int
        Paired queue of reads.

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
    it = iter(readque)
    for x, y in zip(it, it):

        # pre-calculate read metrics
        rid = x & (1 << 30) - 1     # index
        beg = x >> 48               # start coordinate
        end = y >> 48               # end coordinate
        L = (x >> 31) - 1 & 131071  # effective length - 1

        # iterate over genes (ordinal)
        for code in geneque:

            # add gene to cache
            if code & (1 << 31):
                genes[code & (1 << 30) - 1] = code >> 48

            # check overlap while removing gene from cache:
            # min(gene end, read end) - max(gene start, read start) >=
            # effective length - 1
            elif (min(code >> 48, end) -
                  max(genes_pop(code & (1 << 30) - 1), beg)) >= L:
                yield rid, code & (1 << 30) - 1


def match_read_gene_quart(geneque, readque):
    """Associate reads with genes by iterating between read region and nearer
    end of gene queue.

    Parameters
    ----------
    geneque : list of int
        Sorted queue of genes.
    readque : list of int
        Paired queue of reads.

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
    over half of the genome, and they cannot be detected without iterating to
    end of the gene queue.

    This method uses bisection to find the insertion point of read start in the
    pre-sorted gene queue, which has O(logn) time.
    """
    n = len(geneque)  # entire search space
    mid = n // 2      # mid point

    # iterate over paired read starts and ends
    it = iter(readque)
    for x, y in zip(it, it):

        # pre-calculate read metrics
        rid = x & (1 << 30) - 1     # index
        beg = x >> 48               # start coordinate
        end = y >> 48               # end coordinate
        L = (x >> 31) - 1 & 131071  # effective length - 1

        # genes starting within read region
        # will record their coordinates
        within = {}
        within_pop = within.pop

        # read starts in left half of gene queue
        if x < geneque[mid]:

            # genes starting before read region
            # no need to record coordinates as overlap is not possible
            before = set()
            before_add = before.add
            before_remove = before.remove

            # locate read start using bisection
            # Python's `bisect` is more efficient than homebrew code
            i = bisect(geneque, x, hi=mid)

            # iterate over gene positions before read region
            for code in geneque[:i]:

                # add gene to cache when it starts
                if code & (1 << 31):
                    before_add(code & (1 << 30) - 1)

                # drop gene from cache when it ends
                else:
                    before_remove(code & (1 << 30) - 1)

            # iterate over gene positions after read start
            for code in geneque[i:]:

                # stop when exceeding read end
                if code > y:
                    break

                # when gene starts, add it and its coordinate to cache
                elif code & (1 << 31):
                    within[code & (1 << 30) - 1] = code >> 48

                # when gene ends, check overlap and remove it from cache
                else:
                    gid = code & (1 << 30) - 1

                    # most genes should start before read region
                    try:
                        before_remove(gid)

                    # if gene started within read region,
                    # overlap = gene end - gene start
                    except KeyError:
                        if (code >> 48) - within_pop(gid) >= L:
                            yield rid, gid

                    # if gene started before read region,
                    # overlap = gene end - read start
                    else:
                        if (code >> 48) - beg >= L:
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

            # locate read start using bisection
            i = bisect(geneque, x, lo=mid + 1, hi=n)

            # iterate over gene positions within read region
            # for i in range(bisect(geneque, x), n):
            while i < n:
                code = geneque[i]

                # stop when exceeding read end
                if code > y:
                    break

                # when gene starts, add it and its coordinate to cache
                if code & (1 << 31):
                    within[code & (1 << 30) - 1] = code >> 48

                # when gene ends, check overlap and remove it from cache
                # gene end must be <= read end
                # if gene not found in cache (gene started before read region),
                # use read start, otherwise use gene start
                elif (code >> 48) - within_pop(code & (1 << 30) - 1, beg) >= L:
                    yield rid, code & (1 << 30) - 1

                # move to next position
                i += 1

            # iterate over gene positions after read region
            for code in geneque[i:]:

                # add gene to cache when it starts
                if code & (1 << 31):
                    after_add(code & (1 << 30) - 1)

                # when gene ends,
                else:

                    # most remaining genes should start after read region
                    try:
                        after_remove(code & (1 << 30) - 1)

                    # otherwise, the gene spans over entire read region
                    except KeyError:

                        # check overlap while removing gene from cache
                        # gene end must be >= read end
                        # if gene not found in cache, use read start, otherwise
                        # use gene start
                        if end - within_pop(code & (1 << 30) - 1, beg) >= L:
                            yield rid, code & (1 << 30) - 1


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
            gid = idmap_[code & (1 << 30) - 1]
            if prefix:
                gid = nucl + gid
            if code >> 31 & 131071:
                res[gid] = 1 - (code >> 48)
            else:
                res[gid] += code >> 48
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
