#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import logging
import click

from woltka.util import readzip
from woltka.parse import infer_align_format, parse_line_ordinal, assign_parser


@click.command()
@click.option(
    '--input-align', '-i', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='input read alignment')
@click.option(
    '--output-map', '-o', required=True,
    type=click.Path(writable=True, dir_okay=False),
    help='output profile')
@click.option(
    '--genetab', '-g', required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='table of gene coordinates on genome sequences')
@click.option(
    '--input-format', '-f', default='auto',
    type=click.Choice(['auto', 'b6o', 'sam', 'map'], case_sensitive=False),
    help='table of gene coordinates on genome sequences')
@click.option(
    '--threshold', '-t', type=click.FLOAT, default=0.8, show_default=True,
    help=('minimum ratio of overlap length vs. alignment length to qualify '
          'for a match'))
@click.option(
    '--ambiguity', '-a', default='norm', show_default=True,
    type=click.Choice(['uniq', 'all', 'norm'], case_sensitive=False),
    help=('how to treat reads that occur multiple times: "uniq": drop '
          'non-unique reads; "all": count each occurrence once; "norm": '
          'count each occurrence 1/k times (k = number of occurrences)'))
@click.option(
    '--lines', '-n', type=click.INT, default=1000000, show_default=True,
    help=('number of lines per chunck for parsing; higher value improves '
          'time efficiency but consumes more memory'))
@click.option(
    '--outmap', '-m', is_flag=True,
    help='generate a read-to-gene(s) map instead of a gene-to-count profile')
@click.option(
    '--prefix', '-p', is_flag=True,
    help='prefix nucleotide ID to gene ID in the format of "nucl_gene"')
def ordinal(input_align, output_map, genetab, input_format, threshold,
            ambiguity, lines, outmap, prefix):
    """Match reads and genes in an ordinal scale on the genome."""

    # config logger
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s %(message)s')

    # output file
    fo = open(output_map, 'w')

    # whether return read map or gene counts (latter saves compute)
    ismap = outmap or input_format == 'map' or ambiguity != 'all'

    # read gene coordinates
    if genetab:
        logging.info('Gene table parsing started.')
        genetab = read_gene_table(readzip(genetab))

        # sort gene coordinates per nucleotide
        for gene, queue in genetab.items():
            genetab[gene] = sorted(queue, key=lambda x: x[0])
        logging.info('Gene table parsing completed.')

    # gene profile to be constructed
    profile = {}

    # identifiers, alignment lengths and coordinates of reads
    # note: for a read, "id" is the identifier (a free text), "idx" is the
    # index (an incremental integer for all reads of the current chunck)
    rids, lenmap, locmap = [], {}, {}

    # read input read map
    logging.info('Input read alignment parsing started.')

    # line counter
    ii = 0

    # associate reads with genes for current chunck
    def _match_reads_genes():
        nonlocal rids, lenmap, locmap
        readmap = {}
        for nucl, loci in locmap.items():

            # merge and sort coordinates
            # question is to merge an unsorted list into a sorted one
            # Python's built-in "timesort" algorithm is efficient at this
            try:
                queue = sorted(genetab[nucl] + loci, key=lambda x: x[0])

            # it's possible that no gene was annotated on the nucleotide
            except KeyError:
                continue

            # map reads to genes
            res = match_read_gene(queue, lenmap[nucl], threshold, ismap)

            # prefix
            pref = nucl if prefix else None

            # merge into master read map (of current chunck)
            if outmap:
                add_match_to_readmap(readmap, res, pref)

            # merge into master profile (of entire sample)
            else:
                add_match_to_profile(profile, res, ismap, pref)

        # write read map
        if outmap:
            logging.info('Mapped reads: {}.'.format(len(readmap)))
            for i, rid in enumerate(rids):
                try:
                    fo.write('{}\t{}\n'.format(rid, ','.join(sorted(
                        readmap[i]))))
                except KeyError:
                    pass

        # free memory
        lenmap, locmap = {}, {}
        logging.info('Parsed {} lines.'.format(ii))

    # parser for read-to-gene(s) map
    def _parse_map_line(line, *args):
        nonlocal profile, rids
        rid, genes = line.rstrip().split('\t')
        rix = len(rids)
        rids.append(rid)
        for gene in genes.split(','):
            profile.setdefault(gene, []).append(rix)

    # read input alignment
    with readzip(input_align) as f:

        # determine alignment format
        if input_format == 'auto':  # auto-determine
            line = f.readline()
            try:
                input_format = infer_align_format(line)
                logging.info('Alignment format: {}.'.format(input_format))
            except ValueError:
                logging.error('Alignment format cannot be determined.')
                sys.exit(1)
            if input_format == 'map':
                ismap = True
            parser = assign_parser(input_format)
            parse_line_ordinal(line, parser, rids, lenmap, locmap)
            ii += 1
        else:
            parser = assign_parser(input_format)

        # parse read map in chuncks
        for line in f:
            parse_line_ordinal(line, parser, rids, lenmap, locmap)
            ii += 1
            if lines and ii % lines == 0:
                _match_reads_genes()
        _match_reads_genes()
        logging.info('Input read alignment parsing completed.')

    # read map is already written
    if outmap:
        logging.info('Output read-to-gene(s) map written.')
        return

    # convert read maps into counts
    if ismap:
        readmap_to_profile(profile, rids, ambiguity == 'norm')
        logging.info('Multiple-occurrence reads processed.')

    # write profile
    for gene, n in sorted(profile.items()):
        fo.write('{}\t{}\n'.format(gene, n))
    logging.info('Output gene profile written.')


def read_gene_table(f):
    """Read coordinates of genes on genomes.

    Parameters
    ----------
    f : file handle
        Gene table file.

    Returns
    -------
    dict of list of tuple of (int, bool, bool, str)
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
    for line in f:
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
            start, end = sorted([int(x[1]), int(x[2])])
            res[nucl].extend((
                (start, True, True, idx),
                (end,  False, True, idx)))
    return res


def match_read_gene(queue, lens, th=0.8, ismap=False):
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
    ismap : bool, optional
        Return read-to-gene map (True) instead of per-gene counts (False).

    Returns
    -------
    dict
        Per-gene counts or read-to-gene map.

    See Also
    --------
    read_gene_table

    Notes
    -----
    This algorithm is the core of the program. It uses a flattened, sorted
    list to store starting and ending coordinates of both genes and reads.
    Only one round of traversal (O(n)) of this list is needed to accurately
    find all gene-read matches.
    """
    match = {}  # read/gene match
    genes = {}  # current genes
    reads = {}  # current reads

    def _add_to_match(rid, gid):
        if ismap:
            match.setdefault(rid, set()).add(gid)
        else:
            match[gid] = match.get(gid, 0) + 1

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
                        _add_to_match(rid, id_)

                # remove it from current genes
                del(genes[id_])

        # the same for reads
        else:
            if is_start:
                reads[id_] = loc
            else:
                for gid, gloc in genes.items():
                    if loc - max(reads[id_], gloc) + 1 >= lens[id_] * th:
                        _add_to_match(id_, gid)
                del(reads[id_])
    return match


def add_match_to_readmap(readmap, match, nucl=None):
    """Merge current read-gene matches to master read map.

    Parameters
    ----------
    readmap : dict
        Master read map.
    match : dict
        Current read map.
    nucl : str, optional
        Prefix nucleotide Id to gene Ids.
    """
    for rid, genes in match.items():
        if nucl:
            genes = {'{}_{}'.format(nucl, x) for x in genes}
        readmap.setdefault(rid, set()).update(genes)


def add_match_to_profile(profile, match, ismap=True, nucl=None):
    """Merge current read-gene matches to master profile.

    Parameters
    ----------
    profile : dict
        Master gene profile.
    match : dict
        Read-gene matches.
    ismap : bool, optional
        Whether matches are a read-to-gene(s) map or simple counts.
    nucl : str, optional
        Prefix nucleotide Id to gene Ids.

    See Also
    --------
    match_read_gene
    """
    # prefix gene Id with nucleotide Id
    def prefix(gene):
        return '{}_{}'.format(nucl, gene) if nucl else gene

    # read-to-gene(s) map
    if ismap:
        for ridx, genes in match.items():
            for gene in genes:
                profile.setdefault(prefix(gene), []).append(ridx)

    # simple counts
    else:
        for gene, count in match.items():
            gene = prefix(gene)
            profile[gene] = profile.get(gene, 0) + count


def readmap_to_profile(profile, rids, normalize=True):
    """Convert read-to-gene map into gene profile.

    Parameters
    ----------
    profile : dict
        Master gene profile.
    rids : list
        Read identifiers.
    normalize : bool, optional
        For a read occurring k times, normalize by k (True), or drop (False).
    """
    # count occurrences of reads
    freqs = {}
    for read in rids:
        freqs[read] = freqs.get(read, 0) + 1

    # treat multiple occurrences
    todel = []
    for gene, ridxx in profile.items():
        n = 0
        for ridx in ridxx:
            k = freqs[rids[ridx]]
            if k > 1:
                if normalize:
                    n += 1 / k
            else:
                n += 1
        n = round(n)
        if n > 0:
            profile[gene] = n
        else:
            todel.append(gene)

    # delete zero values
    for gene in todel:
        del profile[gene]


if __name__ == '__main__':
    ordinal()
