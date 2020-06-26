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


def plain_mapper(fh, fmt=None, n=1000000):
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
    fmt = fmt or infer_align_format(fh)

    # assign parser for given format
    parser = assign_parser(fmt)

    # query queue and subject(s) queue
    qryque, subque = deque(), deque()

    # pre-load method references
    qry_append, sub_append = qryque.append, subque.append

    # parse alignment file
    this = None  # current query Id
    target = n   # target line number at end of current chunk
    for i, line in enumerate(fh):

        # parse current alignment line
        try:
            query, subject = parser(line)[:2]
        except TypeError:
            continue

        # add subject to subject set of the same query Id
        if query == this:
            subque[-1].add(subject)

        # when query Id changes,
        else:

            # line number has reached target
            if i > target:

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
            sub_append({subject})

            # update current query Id
            this = query

    # final flush
    yield qryque, subque


def parse_align_file(fh, mapper, fmt=None, n=None):
    """Read an alignment file in chunks and yield query-to-subject(s) maps.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    mapper : object
        Module for handling alignments.
    fmt : str, optional
        Format of mapping file.
    n : int, optional
        Number of lines per chunck.

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
    n = n or 1000000  # default chunk size: 1 million
    i = 0        # current line number
    j = n        # target line number at end of current chunk
    last = None  # last query Id

    # save attribute lookup overhead
    parse = mapper.parse
    append = mapper.append
    flush = mapper.flush

    # determine file format based on first line
    if fmt is None:
        try:
            line = next(fh)
        except StopIteration:
            return
        fmt = infer_align_format(line)
        parser = assign_parser(fmt)
        try:
            last = parse(line, parser)
            i = 1
        except TypeError:
            pass
        else:
            append()
    else:
        parser = assign_parser(fmt)

    # parse remaining content
    for line in fh:
        try:
            query = parse(line, parser)
            i += 1
        except TypeError:
            continue

        # flush when query Id changes and chunk size was already reached
        if query != last:
            if i >= j:
                yield flush()
                j = i + n
            last = query
        append()

    # finish last chunk
    yield flush()


class Plain(object):
    """Mapping module for plain sequence alignments.

    Attributes
    ----------
    map : dict
        Cache of query-to-subject(s) map.
    buf : (str, str)
        Buffer of last (query, subject) pair.

    See Also
    --------
    parse_align_file
    ordinal.Ordinal

    Notes
    -----
    The design of this class provides an interface for parsing extremely large
    alignment files. The other class of the same type is `Ordinal`. They both
    provide three instance methods: `parse`, `append` and `flush` which can be
    called from parent process.
    """

    def __init__(self):
        """Initiate mapper.
        """
        self.this = None
        self.qryque, self.subque = deque(), deque()

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
            Query and subject cannot be extracted from line.
        """
        self.buf = parser(line)[:2]
        return self.buf[0]

    def append(self):
        """Append buffered last line to cached read map.
        """
        try:
            query, subject = self.buf
        except (AttributeError, TypeError):
            return
        if query == self.this:
            self.subque[-1].add(subject)
        else:
            self.qryque.append(query)
            self.subque.append({subject})
            self.this = query

    def flush(self):
        """Process, return and clear read map.

        Returns
        -------
        deque of str
            Query queue.
        deque of set of str
            Subject(s) queue.
        """
        try:
            return self.qryque, self.subque
        finally:
            self.__init__()


def infer_align_format(fh):
    """Guess the format of an alignment file based on first line.

    Parameters
    ----------
    fh : file handle
        Input alignment file.

    Returns
    -------
    str
        Alignment file format (map, b6o or sam).

    Raises
    ------
    ValueError
        Alignment file cannot be read.
    ValueError
        Alignment file format cannot be determined.

    See Also
    --------
    parse_b6o_line
    parse_sam_line
    """
    # read first line of file
    try:
        line = next(fh)

    # file is empty or not readable
    except StopIteration:
        raise ValueError('Cannot read alignment file.')

    # move pointer back to beginning of file
    fh.seek(0)

    # SAM file header
    if line.split()[0] in ('@HD', '@PG'):
        return 'sam'

    # tab-delimited row
    row = line.rstrip().split('\t')

    # simple query-to-subject map
    if len(row) == 2:
        return 'map'

    # BLAST standard tabular format
    if len(row) >= 12:
        if all(row[i].isdigit() for i in range(3, 10)):
            return 'b6o'

    # SAM format
    if len(row) >= 11:
        if all(row[i].isdigit() for i in (1, 3, 4)):
            return 'sam'

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
    Only first two columns are considered.
    """
    try:
        query, subject = line.rstrip().split('\t', 2)[:2]
        return query, subject
    except ValueError:
        pass


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
    x = line.rstrip().split('\t')
    qseqid, sseqid, length, score = x[0], x[1], int(x[3]), float(x[11])
    sstart, send = sorted([int(x[8]), int(x[9])])
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
    x = line.rstrip().split('\t')
    qname, rname = x[0], x[2]  # query and subject identifiers

    # skip unmapped
    if rname == '*':
        return
    pos = int(x[3])  # leftmost mapping position

    # parse CIGAR string
    length, offset = cigar_to_lens(x[5])

    # append strand to read Id if not already
    if not qname.endswith(('/1', '/2')):
        flag = int(x[1])
        if flag & (1 << 6):  # forward strand: bit 64
            qname += '/1'
        elif flag & (1 << 7):  # reverse strand: bit 128
            qname += '/2'

    return qname, rname, None, length, pos, pos + offset - 1


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
    x = line.rstrip().split('\t')
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
        return None
    x = line.rstrip().split('\t')
    return x[0], x[1], int(x[3]), int(x[5])
