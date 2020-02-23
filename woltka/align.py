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


def parse_align_file(fh, proc, fmt=None, n=1000000):
    """Read an alignment file in chunks and yield query-to-subject(s) maps.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    proc : object
        Module for handling alignments.
    fmt : str
        Format of mapping file.
    n : int, optional
        Number of lines per chunck.

    Yields
    ------
    dict of list
        Query-to-subject(s) map.

    Notes
    -----
    The design of this function aims to couple with the extremely large size of
    typical alignment files. It reads the entire file sequentially, pauses and
    processes current cache for every _n_ lines, yields and clears cache, then
    proceeds.
    """
    i = 0        # current line number
    j = n        # target line number at end of current chunk
    last = None  # last query Id

    # determine file format based on first line
    if fmt is None:
        try:
            line = next(fh)
        except StopIteration:
            return
        fmt = infer_align_format(line)
        parser = assign_parser(fmt)
        try:
            last = proc.parse(line, parser)
            i = 1
        except TypeError:
            pass
        proc.append()
    else:
        parser = assign_parser(fmt)

    # parse remaining content
    for line in fh:
        try:
            query = proc.parse(line, parser)
            i += 1
        except TypeError:
            continue

        # flush when query Id changes and chunk size was already reached
        if query != last:
            if i >= j:
                yield proc.flush()
                j = i + n
            last = query
        proc.append()

    # finish last chunk
    yield proc.flush()


class Plain(object):
    """Processor of plain sequence alignments.

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
        """Initiate processor.
        """
        self.map = {}

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
        self.map.setdefault(query, []).append(subject)

    def flush(self):
        """Process, return and clear read map.

        Returns
        -------
        Processed read map.
        """
        res, self.map = self.map, {}
        return res


def read_align(fh, fmt=None):
    """Read an alignment file into a query-to-subject(s) map.

    Parameters
    ----------
    fh : file handle
        Mapping file to read.
    fmt : str
        Format of mapping file.

    Returns
    -------
    dict of list
        Query-to-subject(s) map.

    Notes
    -----
    When one query occurs multiple times, all its subjects, including
    duplicates, are added to the result.
    """
    res = {}

    # determine file format based on first line
    if fmt is None:  # auto-determine
        line = fh.readline()
        fmt = infer_align_format(line)
        parser = assign_parser(fmt)
        try:
            query, subject = parser(line)[:2]
        except TypeError:
            pass
        else:
            res.setdefault(query, []).append(subject)
    else:
        parser = assign_parser(fmt)

    # parse remaining content
    for line in fh:
        try:
            query, subject = parser(line)[:2]
        except TypeError:
            continue
        else:
            res.setdefault(query, []).append(subject)

    return res


def infer_align_format(line):
    """Guess the format of an alignment file based on first line.

    Parameters
    ----------
    line : str
        First line of alignment file.

    Returns
    -------
    str or None
        Alignment file format (map, b6o or sam).

    Raises
    ------
    ValueError
        Alignment file format cannot be determined.

    See Also
    --------
    parse_b6o_line
    parse_sam_line
    """
    if line.split()[0] == '@HD':
        return 'sam'
    row = line.rstrip().split('\t')
    if len(row) == 2:
        return 'map'
    if len(row) >= 12:
        if all(row[i].isdigit() for i in range(3, 10)):
            return 'b6o'
    if len(row) >= 11:
        if all(row[i].isdigit() for i in (1, 3, 4)):
            return 'sam'
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
        raise ValueError(f'Invalid format code: {fmt}.')


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
    str, str
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
    """Parse a line in a SAM format (sam).

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
    if line.startswith('@'):
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
    int, int
        Alignment length.
        Offset in subject sequence.

    Raises
    ------
    ValueError
        CIGAR string is missing.
    """
    if cigar in ('', '*'):
        raise ValueError('Missing CIGAR string.')
    align, offset = 0, 0
    ops = 'DHIMNPSX='
    n = ''  # current step size
    for c in cigar:
        if c in ops:
            if c in 'M=X':
                align += int(n)
            if c in 'MDN=X':
                offset += int(n)
            n = ''
        else:
            n += c
    return align, offset


def parse_kraken(line):
    """Parse a line in a Kraken mapping file.

    Parameters
    ----------
    line : str
        Line to parse.

    Returns
    -------
    tuple of (str, str)
        Query, subject.

    Notes
    -----
    Kraken2 output format:
        C/U, sequence Id, taxonomy Id, length, LCA mapping

    .. _Kraken2 manual:
        https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
    """
    x = line.rstrip().split('\t')
    return x[1], x[2] if x[0] == 'C' else None, None


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
    return x[0], x[1], int(x[2]), int(x[3])
