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

    this = None  # current query Id
    target = n   # target line number at end of current chunk

    # parse alignment file
    for i, (query, subject) in enumerate(parser(chain(iter(head), fh))):

        # add subject to subject set of the same query Id
        if query == this:
            subque[-1].add(subject)

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
            sub_append({subject})

            # update current query Id
            this = query

    # final flush
    yield qryque, subque


def range_mapper(fh, fmt=None, n=1000):
    """Read an alignment file and yield maps of query to subject(s) and their
    ranges.

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
    deque of dict of str to list of int
        Subject-to-ranges queue.

    Notes
    -----
    Same as `plain_mapper`, except that it also returns subject ranges.

    Ranges are stored as a one-dimensional, interleaved list of start1, end1,
    start2, end2, start3, end3...

    See Also
    --------
    plain_mapper
    """
    fmt, head = (fmt, []) if fmt else infer_align_format(fh)
    parser = assign_parser(fmt, ext=True)
    qryque, subque = deque(), deque()
    qry_append, sub_append = qryque.append, subque.append
    this = None
    target = n
    for i, row in enumerate(parser(chain(iter(head), fh))):
        query, subject, _, _, start, end = row[:6]

        # range must be positive integers
        if start and end:

            if query == this:
                subque[-1].setdefault(subject, []).extend((start, end))
            else:
                if i >= target:
                    yield qryque, subque
                    qryque, subque = deque(), deque()
                    qry_append, sub_append = qryque.append, subque.append
                    target = i + n
                qry_append(query)

                # return subject Id and range
                sub_append({subject: [start, end]})

                this = query
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
    parse_b6o_file
    parse_sam_file

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


def assign_parser(fmt, ext=False):
    """Assign parser function based on format code.

    Parameters
    ----------
    fmt : str
        Alignment file format code.
    ext : bool, optional
        Whether to get extra information.

    Returns
    -------
    function
        Alignment parser function.
    """
    if fmt == 'map':  # simple map of query <tab> subject
        return parse_map_file
    if fmt == 'b6o':  # BLAST format
        return parse_b6o_file_ext if ext else parse_b6o_file
    elif fmt == 'sam':  # SAM format
        return parse_sam_file_ext if ext else parse_sam_file
    else:
        raise ValueError(f'Invalid format code: "{fmt}".')


def parse_map_file(fh, *args):
    """Parse a simple mapping file.

    Parameters
    ----------
    fh : file handle
        Mapping file to parse.
    args : list, optional
        Placeholder for caller compatibility.

    Yields
    ------
    tuple of (str, str)
        Query and subject.

    Notes
    -----
    Only the first two columns are considered.
    """
    for line in fh:
        query, found, rest = line.partition('\t')
        if found:
            subject, found, rest = rest.partition('\t')
            yield query, subject.rstrip()


def parse_b6o_file(fh):
    """Parse a BLAST tabular file (b6o) to get basic information.

    Parameters
    ----------
    fh : file handle
        BLAST tabular file to parse.

    Yields
    ------
    tuple of (str, str)
        Query, subject.

    Notes
    -----
    BLAST tabular format:
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send
        evalue bitscore

    .. _BLAST manual:
        https://www.ncbi.nlm.nih.gov/books/NBK279684/
    """
    for line in fh:
        try:
            qseqid, sseqid, _ = line.split('\t', 2)
        except ValueError:
            continue
        else:
            yield qseqid, sseqid


def parse_b6o_file_ext(fh):
    """Parse a BLAST tabular file (b6o) to get extra information.

    Parameters
    ----------
    fh : file handle
        BLAST tabular file to parse.

    Yields
    ------
    tuple of (str, str, float, int, int, int)
        Query, subject, score, length, start, end.
    """
    for line in fh:
        x = line.split('\t')
        try:
            qseqid, sseqid, length, score = x[0], x[1], int(x[3]), float(x[11])
        except IndexError:
            continue
        sstart, send = sorted((int(x[8]), int(x[9])))
        yield qseqid, sseqid, score, length, sstart, send


def parse_sam_file(fh):
    """Parse a SAM file (sam) to get basic information.

    Parameters
    ----------
    fh : file handle
        SAM file to parse.

    Yields
    ------
    tuple of (str, str)
        Query, subject.

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

    This is by far the fastest solution. Several other solutions were tested,
    including readlines(chunk), csv.reader, and pd_read_csv(chunk). They were
    slower.
    """
    last = None

    # parse head
    for line in fh:

        # currently, it just skips head
        if line[0] != '@':
            last = line
            break

    # include current line
    it = chain([last], fh) if last else fh

    # parse body
    for line in it:

        # relevant fields
        qname, flag, rname, _ = line.split('\t', 3)

        # skip unmapped
        if rname == '*':
            continue

        # append strand to read Id if not already
        if qname[-2:] not in ('/1', '/2'):
            flag = int(flag)

            # forward strand: bit 64
            if flag & (1 << 6):
                qname += '/1'

            # reverse strand: bit 128
            elif flag & (1 << 7):
                qname += '/2'

        yield qname, rname


def parse_sam_file_ext(fh):
    """Parse a SAM file (sam) to get extra information.

    Parameters
    ----------
    fh : file handle
        SAM file to parse.

    Yields
    ------
    tuple of (str, str, None, int, int, int)
        Query, subject, None, length, start, end.
    """
    last = None
    for line in fh:
        if line[0] != '@':
            last = line
            break

    it = chain([last], fh) if last else fh
    for line in it:

        # relevant fields
        qname, flag, rname, pos, _, cigar, _ = line.split('\t', 6)
        if rname == '*':
            continue

        # leftmost mapping position
        pos = int(pos)

        # parse CIGAR string
        length, offset = cigar_to_lens(cigar)

        if qname[-2:] not in ('/1', '/2'):
            flag = int(flag)
            if flag & (1 << 6):
                qname += '/1'
            elif flag & (1 << 7):
                qname += '/2'

        yield qname, rname, None, length, pos, pos + offset - 1


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


def cigar_to_lens_ord(cigar):
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
    This is an alternative solution based on unicode numbers of characters.
    It is slower than the currently adopted solution. But since it is purely
    based on numbers, there may be room for future optimization.
    """
    align, offset = 0, 0
    n = 0
    for c in map(ord, cigar):
        if c < 60:
            n = n * 10 + c - 48
        else:
            if c in (77, 88, 61):
                align += n
            elif c in (68, 78):
                offset += n
            n = 0
    return align, align + offset


def parse_kraken(fh):
    """Parse a Kraken mapping file.

    Parameters
    ----------
    fh : file handle
        Kraken mapping file to parse.

    Yields
    ------
    tuple of (str, str)
        Query and subject.

    Notes
    -----
    Kraken2 output format:
        C/U, sequence ID, taxonomy ID, length, LCA mapping

    .. _Kraken2 manual:
        https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
    """
    for line in fh:
        x = line.split('\t')
        if x[0] == 'C':
            yield x[1], x[2]


def parse_centrifuge(fh):
    """Parse a Centrifuge mapping file.

    Parameters
    ----------
    fh : file handle
        Centrifuge mapping file to parse.

    Yields
    ------
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
    for line in fh:
        if not line.startswith('readID'):
            x = line.split('\t')
            yield x[0], x[1], int(x[3]), int(x[5])


def parse_sam_file_pd(fh, n=65536):
    """Parse a SAM file (sam) using Pandas.

    Parameters
    ----------
    fh : file handle
        SAM file to parse.
    n : int, optional
        Chunk size.

    Yields
    ------
    None

    Notes
    -----
    This is a SAM file parser using Pandas. It is slower than the current
    parser. The `read_csv` is fast, but the data frame manipulation slows
    down the process. It is here for reference only.
    """
    return
    # with pd.read_csv(fp, sep='\t',
    #                  header=None,
    #                  comment='@',
    #                  na_values='*',
    #                  usecols=[0, 1, 2, 3, 5],
    #                  names=['qname', 'flag', 'rname', 'pos', 'cigar'],
    #                  dtype={'qname': str,
    #                         'flag':  np.uint16,
    #                         'rname': str,
    #                         'pos':   int,
    #                         'cigar': str},
    #                  chunksize=n) as reader:
    #     for chunk in reader:
    #         chunk.dropna(subset=['rname'], inplace=True)
    #         # this is slow, because of function all
    #         chunk['length'], offset = zip(*chunk['cigar'].apply(
    #             cigar_to_lens))
    #         chunk['right'] = chunk['pos'] + offset - 1
    #         # this is slow, because of function all
    #         # chunk['qname'] = chunk[['qname', 'flag']].apply(
    #         #   qname_by_flag, axis=1)
    #         # this is a faster method
    #         chunk['qname'] += np.where(
    #             chunk['qname'].str[-2:].isin(('/1', '/2')), '',
    #             np.where(np.bitwise_and(chunk['flag'], (1 << 6)), '/1',
    #                      np.where(np.bitwise_and(chunk['flag'], (1 << 7)),
    #                      '/2', '')))
    #         chunk['score'] = 0
    #         yield from chunk[['qname', 'rname', 'score', 'length',
    #                           'pos', 'right']].values
