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

from itertools import chain


class Ordinal(object):
    """Processor for matching reads and genes in an ordinal system.

    Attributes
    ----------
    coords : dict
        Gene coordinates table.
    prefix : bool
        Prefix gene IDs with Nucleotide IDs.
    th : float
        Minimum threshold of overlap length : alignment length for a match.
    buf : list
        Buffer of last alignment.
    rids : list
        Read identifiers in cache.
    lenmap : dict
        Alignment length map in cache.
    locmap : dict
        Alignment location map in cache.

    See Also
    --------
    align.Plain

    Notes
    -----
    For a read, "id" is the identifier (a free text), "idx" is the index (an
    incremental integer for all reads of the current chunk).
    """

    def __init__(self, coords, prefix=False, th=0.8):
        """Initiate processor.

        Parameters
        ----------
        coords : dict
            Gene coordinates table.
        prefix : bool
            Whether prefix gene IDs.
        th : float
            Overlap : alignment threshold.
        """
        self.coords = coords
        self.prefix = prefix
        self.th = th
        self.buf = None
        self.clear()

    def parse(self, line, parser):
        """Parse one line in alignment file.

        Parameters
        ----------
        line : str
            Line in alignment file.
        parser : callable
            Function to parse the line.

        Returns
        -------
        str
            Query identifier.

        Raises
        ------
        TypeError
            Necessary information cannot be extracted from line.
        """
        self.buf = parser(line)[:6]
        return self.buf[0]

    def append(self):
        """Append buffered last line to cached alignment map.

        See Also
        --------
        read_gene_coords
        """
        try:
            query, subject, _, length, start, end = self.buf
        except (AttributeError, TypeError):
            return
        idx = len(self.rids)
        self.rids.append(query)
        self.lenmap.setdefault(subject, {})[idx] = length
        self.locmap.setdefault(subject, []).extend((
            (start, True, False, idx),
            (end,  False, False, idx)))

    def flush(self):
        """Process, return and clear read map.

        Returns
        -------
        deque of str
            Query queue.
        deque of set of str
            Subject(s) queue.
        """
        from_iterable = chain.from_iterable
        th = self.th
        rids = self.rids
        coords = self.coords
        lenmap = self.lenmap
        prefix = self.prefix

        matches = []
        match_append = matches.append
        for nucl, loci in self.locmap.items():

            # merge and sort coordinates
            # question is to merge an unsorted list into a sorted one
            # Python's built-in "timesort" algorithm is efficient at this
            try:
                queue = sorted(coords[nucl] + loci, key=lambda x: x[0])

            # it's possible that no gene was annotated on the nucleotide
            except KeyError:
                continue

            # map reads to genes
            match = match_read_gene(queue, lenmap[nucl], th)

            # add nucleotide prefix
            if prefix:
                for idx, genes in match.items():
                    match[idx] = [f'{nucl}_{x}' for x in genes]

            # add to matches
            match_append(match)

        # generate query and subject(s) queues
        # qryque, subque = deque(), deque()
        # qry_append, sub_append = qryque.append, subque.append

        # for idx in set(from_iterable(matches)):
        #     qry_append(rids[idx])
        #     sub_append(set(from_iterable(x.get(idx, []) for x in matches)))

        # merge and de-duplicate matches
        idxdic = {idx: i for i, idx in enumerate(sorted(set(from_iterable(
            matches))))}
        subque = [set() for _ in idxdic]
        for match in matches:
            for idx, genes in match.items():
                subque[idxdic[idx]].update(genes)
        qryque = [rids[idx] for idx in idxdic]

        try:
            return qryque, subque
        finally:
            self.clear()

    def clear(self):
        self.rids = []
        self.lenmap = {}
        self.locmap = {}


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
    map_read_gene

    Notes
    -----
    This data structure is central to this algorithm. Starting and ending
    coordinates of each gene are separated and flattened into a sorted list.
    which enables only one round of list traversal for the entire set of genes
    plus reads.
    """
    res = {}
    nucl = None
    for line in fh:
        line = line.rstrip()

        # ">" or "#" indicates nucleotide name
        if line.startswith(('>', '#')):
            nucl = line[1:].strip()

            # double ">" or "#" indicates genome name
            if not nucl.startswith(('>', '#')):
                res[nucl] = []
        else:
            x = line.split('\t')
            idx = x[0]

            # start and end are based on genome, not gene itself
            try:
                start, end = sorted([int(x[1]), int(x[2])])
            except (IndexError, ValueError):
                raise ValueError(
                    f'Cannot extract coordinates from line: "{line}".')
            res[nucl].extend((
                (start, True, True, idx),
                (end,  False, True, idx)))

    # sort gene coordinates per nucleotide
    if sort:
        for nucl, queue in res.items():
            res[nucl] = sorted(queue, key=lambda x: x[0])

    return res


def whether_prefix(coords):
    """determine whether gene IDs should be prefixed with nucleotide IDs.

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
        for coord, is_start, is_gene, id_ in queue:
            try:
                if genes[id_] == is_start:
                    return True
            except KeyError:
                genes[id_] = is_start
    return False


def match_read_gene(queue, lens, th=0.8):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : list of tuple
        Sorted list of elements.
        (loc, is_start, is_gene, id)
    lens : dict
        Read-to-alignment length map.
    th : float, optional
        Threshold for read/gene overlapping fraction.

    Returns
    -------
    dict
        Read-to-gene(s) map.

    See Also
    --------
    read_gene_coords

    Notes
    -----
    This algorithm is the core of this module. It uses a flattened, sorted
    list to store starting and ending coordinates of both genes and reads.
    Only one round of traversal (O(n)) of this list is needed to accurately
    find all gene-read matches.
    """
    match = {}  # read/gene match
    genes = {}  # current genes
    reads = {}  # current reads

    # cache function reference
    add_to_match = match.setdefault

    # walk through flattened queue of reads and genes
    for loc, is_start, is_gene, id_ in queue:
        if is_gene:

            # when a gene starts, added to current genes
            if is_start:
                genes[id_] = loc

            # when a gene ends,
            else:

                # check current reads
                for rid, rloc in reads.items():

                    # add to match if read/gene overlap is long enough
                    if loc - max(genes[id_], rloc) + 1 >= lens[rid] * th:
                        add_to_match(rid, []).append(id_)

                # remove it from current genes
                del(genes[id_])

        # the same for reads
        else:
            if is_start:
                reads[id_] = loc
            else:
                for gid, gloc in genes.items():
                    if loc - max(reads[id_], gloc) + 1 >= lens[id_] * th:
                        add_to_match(id_, []).append(gid)
                del(reads[id_])
    return match
