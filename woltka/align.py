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

from .util import is_int, is_float, is_pos_int


def plain_mapper(fh, fmt=None, n=1000):
    """Read an alignment file in chunks and yield query-to-subject(s) maps.

    Parameters
    ----------
    fh : file handle
        Alignment file to parse.
    fmt : str, optional
        Alignment file format.
    n : int, optional
        Number of unique queries per chunk.

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
    processes current cache for every _n_ unique queries, yields and clears
    cache, then proceeds.
    """
    # determine alignment file format
    if not fmt:
        fmt, head = infer_align_format(fh)
        it = chain(iter(head), fh)
    else:
        it = fh

    # assign parser for given format
    parser = assign_parser(fmt)

    # query queue and subject(s) queue
    qryque, subque = deque(), deque()

    # pre-load method references
    qryapd, subapd = qryque.append, subque.append

    # target query count at end of current chunk
    target = n

    # parse alignment file
    for i, (query, subjects) in enumerate(parser(it)):

        # line number has reached target
        if i == target:

            # flush current queues of queries and subjects
            yield qryque, subque

            # reset queues and method references
            qryque, subque = deque(), deque()
            qryapd, subapd = qryque.append, subque.append

            # next target query count
            target += n

        # append query and subjects
        qryapd(query)
        subapd(set(subjects))

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
        Number of unique queries per chunk.

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
    .coverage.merge_ranges
    """
    if not fmt:
        fmt, head = infer_align_format(fh)
        it = chain(iter(head), fh)
    else:
        it = fh
    parser = assign_parser(fmt, ext=True)
    qryque, subque = deque(), deque()
    qryapd, subapd = qryque.append, subque.append

    target = n
    for i, (query, records) in enumerate(parser(it)):
        if i == target:
            yield qryque, subque
            qryque, subque = deque(), deque()
            qryapd, subapd = qryque.append, subque.append
            target += n

        # generate a mapping of subjects to interleaved starts and ends
        ranges = {}
        for subject, _, _, start, end in records:

            # start and end must be positive integers
            if start and end:

                # combine ranges on the same subject
                ranges.setdefault(subject, []).extend((start, end))

        # append query and ranges
        if ranges:
            qryapd(query)
            subapd(ranges)

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
        Alignment file format (map, b6o, sam or paf).
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
    parse_map_file
    parse_b6o_file
    parse_sam_file
    parse_paf_file

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

    if len(row) >= 12:

        # BLAST standard tabular format
        if all(row[i].isdigit() for i in range(3, 10)):
            return 'b6o', [line]

        # PAF tabular format
        elif row[4] in '+-' and all(row[i].isdigit() for i in
                                    (1, 2, 3, 6, 7, 8, 9, 10, 11)):
            return 'paf', [line]

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
    callable
        Alignment parser function.
    """
    if fmt == 'map':  # simple map of query <tab> subject
        return parse_map_file
    if fmt == 'b6o':  # BLAST format
        return parse_b6o_file_ext if ext else parse_b6o_file
    elif fmt == 'sam':  # SAM format
        return parse_sam_file_ext if ext else parse_sam_file
    elif fmt == 'paf':  # PAF format
        return parse_paf_file_ext if ext else parse_paf_file
    else:
        raise ValueError(f'Invalid format code: "{fmt}".')


def parse_sam_file(fh):
    """Parse a SAM file (sam) to get basic information.

    Parameters
    ----------
    fh : file handle
        SAM file to parse.

    Yields
    ------
    str
        Query.
    list of str
        Subjects.

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
        https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output

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

    # include last line in body
    it = chain([last], fh) if last else fh

    # current qname
    this = None

    # current subjects by mate: unpaired, forward, reverse
    pool = ([], [], [])

    # parse body
    for line in it:

        # relevant fields
        qname, flag, rname, _ = line.split('\t', 3)

        # skip unmapped
        if rname == '*':
            continue

        # extract mate from flag (forward: bit 6; reverse: bit 7)
        mate = int(flag) >> 6 & 3

        # previous query completes
        if qname != this:

            # yield query and subjects
            if pool[0]:
                yield this, pool[0]
            if pool[1]:
                yield this + '/1', pool[1]
            if pool[2]:
                yield this + '/2', pool[2]

            # reset query and subjects
            this, pool = qname, ([], [], [])

        # add subject to pool
        pool[mate].append(rname)

    # final yield
    if pool[0]:
        yield this, pool[0]
    if pool[1]:
        yield this + '/1', pool[1]
    if pool[2]:
        yield this + '/2', pool[2]


def parse_sam_file_ext(fh):
    """Parse a SAM file (sam) to get extra information.

    Parameters
    ----------
    fh : file handle
        SAM file to parse.

    Yields
    ------
    str
        Query.
    list of (str, None, int, int, int)
        Records of (subject, None, length, start, end).
    """
    last = None
    for line in fh:
        if line[0] != '@':
            last = line
            break

    it = chain([last], fh) if last else fh

    this, pool = None, ([], [], [])
    for line in it:

        # relevant fields
        qname, flag, rname, pos, _, cigar, _ = line.split('\t', 6)
        if rname == '*':
            continue
        mate = int(flag) >> 6 & 3

        # leftmost mapping position
        pos = int(pos)

        # parse CIGAR string
        length, offset = cigar_to_lens(cigar)

        # yield and reset
        if qname != this:
            if pool[0]:
                yield this, pool[0]
            if pool[1]:
                yield this + '/1', pool[1]
            if pool[2]:
                yield this + '/2', pool[2]
            this, pool = qname, ([], [], [])

        # append
        pool[mate].append((rname, None, length, pos, pos + offset - 1))

    # final yield
    if pool[0]:
        yield this, pool[0]
    if pool[1]:
        yield this + '/1', pool[1]
    if pool[2]:
        yield this + '/2', pool[2]


def check_sam_file(fh):
    """Check if a SAM file (sam) is valid.

    Parameters
    ----------
    fh : file handle
        SAM file to check.

    Raises
    ------
    ValueError
        If the file format is invalid.

    Notes
    -----
    This function only checks if the file can be parsed by this program. It
    does not perform a full-scale check following the format's definition.
    """
    count = 0
    for line in fh:
        if line.startswith('@'):
            continue
        line = line.rstrip()
        x = line.split('\t')
        if not (len(x) >= 11
                and x[0]              # QNAME: non-empty string
                and is_int(x[1])      # FLAG: non-negative integer
                and x[2]              # RNAME: non-empty string
                and is_pos_int(x[3])  # POS: positive integer
                and is_cigar(x[5])    # CIGAR: valid CIGAR string
                ):
            raise ValueError(f'Invalid SAM alignment: "{line}".')
        count += 1
    if count == 0:
        raise ValueError('No SAM alignment is found in the file.')


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


def is_cigar(cigar):
    """Check if a CIGAR string is valid.

    Parameters
    ----------
    cigar : str
        CIGAR string to check.

    Returns
    -------
    bool
        Whether the CIGAR string is valid.
    """
    if cigar == '*':
        return True
    is_code = None
    for c in cigar:
        if c in '0123456789':
            if is_code is not False:
                is_code = False
        elif c in 'MDIHNPSX=':
            if is_code is not False:
                return False
            is_code = True
        else:
            return False
    if is_code is not True:
        return False
    return True


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
    str
        Query.
    list of str
        Subjects.

    Notes
    -----
    Only the first two columns are considered.
    """
    # current query and subjects
    this, pool = None, []

    # find first query (micro-optimization)
    for line in fh:

        # extract 1st and 2nd fields (query and subject)
        query, found, rest = line.partition('\t')
        if found:
            subject, _, _ = rest.partition('\t')

            # initialize query
            this, pool = query, [subject.rstrip()]
            break

    # iterate over the rest of file
    for line in fh:
        query, found, rest = line.partition('\t')
        if found:
            subject, _, _ = rest.partition('\t')

            # if query is the same, add subject to pool
            if query == this:
                pool.append(subject.rstrip())

            # if not, yield and reset query and subjects
            else:
                yield this, pool
                this, pool = query, [subject.rstrip()]

    # final yield
    if this is not None:
        yield this, pool


def check_map_file(fh):
    """Check if simple mapping file is valid.

    Parameters
    ----------
    fh : file handle
        Mapping file to check.

    Raises
    ------
    ValueError
        If the file format is invalid.
    """
    count = 0
    for line in fh:
        if line.startswith('#'):
            continue
        x = line.rstrip().split('\t')
        if len(x) >= 2 and x[0] and x[1]:
            count += 1
    if count == 0:
        raise ValueError('No mapping is found in the file.')


def parse_b6o_file(fh):
    """Parse a BLAST tabular file (b6o) to get basic information.

    Parameters
    ----------
    fh : file handle
        BLAST tabular file to parse.

    Yields
    ------
    str
        Query.
    list of str
        Subjects.

    Notes
    -----
    BLAST tabular format:
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send
        evalue bitscore

    .. _BLAST manual:
        https://www.ncbi.nlm.nih.gov/books/NBK279684/
    """
    this, pool = None, []

    # head
    for line in fh:
        try:
            qseqid, sseqid, _ = line.split('\t', 2)
        except ValueError:
            continue
        this, pool = qseqid, [sseqid]
        break

    # body
    for line in fh:
        try:
            qseqid, sseqid, _ = line.split('\t', 2)
        except ValueError:
            continue
        if qseqid == this:
            pool.append(sseqid)
        else:
            yield this, pool
            this, pool = qseqid, [sseqid]

    # foot
    if this is not None:
        yield this, pool


def parse_b6o_file_ext(fh):
    """Parse a BLAST tabular file (b6o) to get extra information.

    Parameters
    ----------
    fh : file handle
        BLAST tabular file to parse.

    Yields
    ------
    str
        Query.
    list of (str, float, int, int, int)
        Records of (subject, score, length, start, end).
    """
    this, pool = None, []

    # head
    for line in fh:
        x = line.split('\t')
        try:
            qseqid, sseqid, length, score = x[0], x[1], int(x[3]), float(x[11])
        except IndexError:
            continue
        sstart, send = sorted((int(x[8]), int(x[9])))
        this, pool = qseqid, [(sseqid, score, length, sstart, send)]
        break

    # body
    for line in fh:
        x = line.split('\t')
        try:
            qseqid, sseqid, length, score = x[0], x[1], int(x[3]), float(x[11])
        except IndexError:
            continue
        sstart, send = sorted((int(x[8]), int(x[9])))
        if qseqid == this:
            pool.append((sseqid, score, length, sstart, send))
        else:
            yield this, pool
            this, pool = qseqid, [(sseqid, score, length, sstart, send)]

    # foot
    if this is not None:
        yield this, pool


def check_b6o_file(fh):
    """Check if a BLAST tabular file (b6o) is valid.

    Parameters
    ----------
    fh : file handle
        BLAST tabular file to check.

    Raises
    ------
    ValueError
        If the file format is invalid.

    Notes
    -----
    This function only checks if the file can be parsed by this program. It
    does not perform a full-scale check following the format's definition.
    """
    count = 0
    for line in fh:
        if line.startswith('#'):
            continue
        line = line.rstrip()
        x = line.split('\t')
        if not (len(x) >= 12
                and x[0]              # qseqid: non-empty string
                and x[1]              # sseqid: non-empty string
                and is_pos_int(x[3])  # length: positive integer
                and is_int(x[8])      # sstart: integer
                and is_int(x[9])      # send: integer
                and is_float(x[11])   # bitscore: float
                ):
            raise ValueError(f'Invalid BLAST hit record: "{line}".')
        count += 1
    if count == 0:
        raise ValueError('No BLAST hit record is found in the file.')


def parse_paf_file(fh):
    """Parse a PAF file (paf) to get basic information.

    Parameters
    ----------
    fh : file handle
        PAF file to parse.

    Yields
    ------
    str
        Query.
    list of str
        Subjects.

    Notes
    -----
    PAF format (first 12 columns):
        1 	string 	Query sequence name
        2 	int 	Query sequence length
        3 	int 	Query start (0-based; BED-like; closed)
        4 	int 	Query end (0-based; BED-like; open)
        5 	char 	Relative strand: "+" or "-"
        6 	string 	Target sequence name
        7 	int 	Target sequence length
        8 	int 	Target start on original strand (0-based)
        9 	int 	Target end on original strand (0-based)
        10 	int 	Number of residue matches
        11 	int 	Alignment block length
        12 	int 	Mapping quality (0-255; 255 for missing)

    .. _PAF format specification:
        https://github.com/lh3/miniasm/blob/master/PAF.md
    """
    this, pool = None, []

    # head
    for line in fh:
        try:
            qname, _, _, _, _, tname, _ = line.split('\t', 6)
        except ValueError:
            continue
        this, pool = qname, [tname]
        break

    # body
    for line in fh:
        try:
            qname, _, _, _, _, tname, _ = line.split('\t', 6)
        except ValueError:
            continue
        if qname == this:
            pool.append(tname)
        else:
            yield this, pool
            this, pool = qname, [tname]

    # foot
    if this is not None:
        yield this, pool


def parse_paf_file_ext(fh):
    """Parse a PAF file (paf) to get extra information.

    Parameters
    ----------
    fh : file handle
        PAF file to parse.

    Yields
    ------
    str
        Query.
    list of (str, int, int, int, int)
        Records of (subject, score, length, start, end).
    """
    this, pool = None, []

    # head
    for line in fh:
        x = line.split('\t')
        try:
            rec = (x[5], int(x[11]), int(x[10]), int(x[7]) + 1, int(x[8]))
        except (IndexError, ValueError):
            continue
        this, pool = x[0], [rec]
        break

    # body
    for line in fh:
        x = line.split('\t')
        try:
            rec = (x[5], int(x[11]), int(x[10]), int(x[7]) + 1, int(x[8]))
        except (IndexError, ValueError):
            continue
        if x[0] == this:
            pool.append(rec)
        else:
            yield this, pool
            this, pool = x[0], [rec]

    # foot
    if this is not None:
        yield this, pool


def check_paf_file(fh):
    """Check if a PAF file (paf) is valid.

    Parameters
    ----------
    fh : file handle
        PAF file to check.

    Raises
    ------
    ValueError
        If the file format is invalid.

    Notes
    -----
    This function only checks if the file can be parsed by this program. It
    does not perform a full-scale check following the format's definition.
    """
    count = 0
    for line in fh:
        if line.startswith('#'):
            continue
        line = line.rstrip()
        x = line.split('\t')
        if not (len(x) >= 12
                and x[0]               # Query sequence: non-empty string
                and x[5]               # Target sequence: non-empty string
                and is_int(x[7])       # Target start: integer
                and is_int(x[8])       # Target end: integer
                and int(x[8]) > int(x[7])
                and is_pos_int(x[10])  # Alignment length: positive integer
                and is_int(x[11])      # Mapping quality: integer
                ):
            raise ValueError(f'Invalid PAF alignment: "{line}".')
        count += 1
    if count == 0:
        raise ValueError('No PAF alignment is found in the file.')
