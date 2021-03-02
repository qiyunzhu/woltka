#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


"""Functions for parsing alignment / mapping files.

Notes
-----
A parser function operates on a single line in the input file.
A parser function returns a tuple of:
    str : query Id
    str : subject Id
    float : alignment score, optional
    int : alignment length, optional
    start : alignment start (5') coordinate, optional
    start : alignment end (3') coordinate, optional
"""

from collections import deque
from itertools import chain
from functools import lru_cache


def plain_mapper(fh, fmt=None, n=1000):
    """Read an alignment file in chunks and yield query-to-subject(s) maps.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    fmt : str, optional
        Alignment file format.
    n : int, optional
        Number of lines per chunk.

    Yields
    ------
    deque of str
        Query queue.
    deque of set of str
        Subject(s) queue.

    Notes
    -----
    The design of this function aims to couple with the extremely large size of
    typical alignment files. It reads the entire file sequentially, pauses and
    processes current cache for every _n_ lines, yields and clears cache, then
    proceeds.
    """
    # determine alignment file format
    fmt, head = (fmt, []) if fmt else infer_align_format(fh)

    # assign parser for given format
    parser = assign_parser(fmt)

    # query queue and subject(s) queue
    qryque, subque = deque(), deque()

    # pre-load method references
    qry_append, sub_append = qryque.append, subque.append

    # parse alignment file
    this = None  # current query Id
    target = n   # target line number at end of current chunk
    for i, line in enumerate(chain(iter(head), fh)):

        # parse current alignment line
        try:
            query, subject = parser(line)[:2]
        except (TypeError, IndexError):
            continue

        # add subject to subject set of the same query Id
        if query == this:
            subque[-1].append(subject)

        # when query Id changes,
        else:

            # line number has reached target
            if i >= target:

                # flush current queries and subject sets
                # this happens only when: 1) query Id changes, and 2) line
                # number has reached target, so that subjects of the same query
                # won't be separated in multiple flushes
                yield qryque, subque

                # re-initiate queues and method references
                qryque, subque = deque(), deque()
                qry_append, sub_append = qryque.append, subque.append

                # next target line number
                target = i + n

            # create new query and subject set pair
            qry_append(query)
            sub_append([subject])

            # update current query Id
            this = query

    # final flush
    yield qryque, subque


def infer_align_format(fh):
    """Guess the format of an alignment file based on content.

    Parameters
    ----------
    fh : file handle
        Input alignment file.

    Returns
    -------
    str
        Alignment file format (map, b6o or sam).
    list of str
        Lines that are read in order to infer format.

    Raises
    ------
    ValueError
        Alignment file is empty or unreadable.
    ValueError
        Alignment file format cannot be determined.

    See Also
    --------
    parse_b6o_line
    parse_sam_line

    TODO
    ----
    Currently this function guesses format only based on the first line.
    This can be more robust.
    """
    # read first line of file
    try:
        line = next(fh)

    # file is empty or not readable
    except StopIteration:
        raise ValueError('Alignment file is empty or unreadable.')

    # SAM file header
    if line.split()[0] in ('@HD', '@PG'):
        return 'sam', [line]

    # tab-delimited row
    row = line.rstrip().split('\t')

    # simple query-to-subject map
    if len(row) == 2:
        return 'map', [line]

    # BLAST standard tabular format
    if len(row) >= 12:
        if all(row[i].isdigit() for i in range(3, 10)):
            return 'b6o', [line]

    # SAM format
    if len(row) >= 11:
        if all(row[i].isdigit() for i in (1, 3, 4)):
            return 'sam', [line]

    # cannot determine
    raise ValueError('Cannot determine alignment file format.')


def assign_parser(fmt):
    """Assign parser function based on format code.

    Parameters
    ----------
    fmt : str
        Alignment file format code.

    Returns
    -------
    function
        Alignment parser function.
    """
    if fmt == 'map':  # simple map of query <tab> subject
        return parse_map_line
    if fmt == 'b6o':  # BLAST format
        return parse_b6o_line
    elif fmt == 'sam':  # SAM format
        return parse_sam_line
    else:
        raise ValueError(f'Invalid format code: "{fmt}".')


def parse_map_line(line, *args):
    """Parse a line in a simple mapping file.

    Parameters
    ----------
    line : str
        Line to parse.
    args : list, optional
        Placeholder for caller compatibility.

    Returns
    -------
    tuple of (str, str)
        Query and subject.

    Notes
    -----
    Only the first two columns are considered.
    """
    query, found, rest = line.partition('\t')
    if found:
        subject, found, rest = rest.partition('\t')
        return query, subject.rstrip()


def parse_b6o_line(line):
    """Parse a line in a BLAST tabular file (b6o).

    Parameters
    ----------
    line : str
        Line to parse.

    Returns
    -------
    tuple of (str, str, float, int, int, int)
        Query, subject, score, length, start, end.

    Notes
    -----
    BLAST tabular format:
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send
        evalue bitscore

    .. _BLAST manual:
        https://www.ncbi.nlm.nih.gov/books/NBK279684/
    """
    x = line.split('\t')
    qseqid, sseqid, length, score = x[0], x[1], int(x[3]), float(x[11])
    sstart, send = sorted((int(x[8]), int(x[9])))
    return qseqid, sseqid, score, length, sstart, send


def parse_sam_line(line):
    """Parse a line in a SAM file (sam).

    Parameters
    ----------
    line : str
        Line to parse.

    Returns
    -------
    tuple of (str, str, None, int, int, int)
        Query, subject, None, length, start, end.

    Notes
    -----
    SAM format:
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL,
        TAGS

    .. _Wikipedia:
        https://en.wikipedia.org/wiki/SAM_(file_format)
    .. _SAM format specification:
        https://samtools.github.io/hts-specs/SAMv1.pdf
    .. _Bowtie2 manual:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
    """
    # skip header
    if line[0] == '@':
        return

    # relevant fields
    qname, flag, rname, pos, _, cigar, _ = line.split('\t', 6)

    # skip unmapped
    if rname == '*':
        return

    # leftmost mapping position
    pos = int(pos)

    # parse CIGAR string
    length, offset = cigar_to_lens(cigar)

    # append strand to read Id if not already
    if qname[-2:] not in ('/1', '/2'):
        flag = int(flag)

        # forward strand: bit 64
        if flag & (1 << 6):
            qname += '/1'

        # reverse strand: bit 128
        elif flag & (1 << 7):
            qname += '/2'

    return qname, rname, None, length, pos, pos + offset - 1


@lru_cache(maxsize=128)
def cigar_to_lens(cigar):
    """Extract lengths from a CIGAR string.

    Parameters
    ----------
    cigar : str
        CIGAR string.

    Returns
    -------
    int
        Alignment length.
    int
        Offset in subject sequence.

    Notes
    -----
    This function significantly benefits from LRU cache because high-frequency
    CIGAR strings (e.g., "150M") are common and redundant calculations can be
    saved.
    """
    align, offset = 0, 0
    n = ''  # current step size
    for c in cigar:
        if c in 'MDIHNPSX=':
            if c in 'M=X':
                align += int(n)
            elif c in 'DN':
                offset += int(n)
            n = ''
        else:
            n += c
    return align, align + offset


def parse_kraken(line):
    """Parse a line in a Kraken mapping file.

    Parameters
    ----------
    line : str
        Line to parse.

    Returns
    -------
    tuple of (str, str)
        Query and subject.

    Notes
    -----
    Kraken2 output format:
        C/U, sequence ID, taxonomy ID, length, LCA mapping

    .. _Kraken2 manual:
        https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
    """
    x = line.split('\t')
    return (x[1], x[2]) if x[0] == 'C' else (None, None)


def parse_centrifuge(line):
    """Parse a line in a Centrifuge mapping file.

    Parameters
    ----------
    line : str
        Line to parse.

    Returns
    -------
    tuple of (str, str, int, int)
        Query, subject, score, length.

    Notes
    -----
    Centrifuge output format:
        readID, seqID, taxID, score, 2ndBestScore, hitLength, queryLength,
        numMatches

    .. _Centrifuge manual:
        https://ccb.jhu.edu/software/centrifuge/manual.shtml
    """
    if line.startswith('readID'):
        return
    x = line.split('\t')
    return x[0], x[1], int(x[3]), int(x[5])
