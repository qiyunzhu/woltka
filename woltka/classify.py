#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Main module for classification of sequencing data.
"""

from os import makedirs
from os.path import join, basename
import click

from .util import readzip, count_list, add_dict, intize, delnone, allkeys
from .sample import read_ids, match_sample_file, demultiplex
from .align import read_align_file
from .ordinal import read_gene_coords, whether_prefix, ordinal
from .tree import build_tree, find_rank, find_lca


def classify(input_path:  str,
             output_path: str,
             input_fmt:   str = None,
             input_ext:   str = None,
             sample_ids: list = None,
             demux:      bool = None,
             rank_lst:    str = None,
             multi:       str = True,
             lca:        bool = True,
             ixend:      bool = False,
             coords_fp:   str = None,
             map_fps:    list = None,
             names_fp:    str = None,
             nodes_fp:    str = None,
             newick_fp:   str = None,
             ranktb_fp:   str = None,
             lineage_fp:  str = None) -> dict:
    """Main classification workflow.

    Parameters
    ----------
    input_path : str
        Path to input alignment file or directory of alignment files.
            If file,

            If directory, demultiplexing is off unless specified.
    output_path : str
        Path to output profile file or directory of profile files.
    input_fmt : str, optional
        Format of input alignment file. Options:
            'b6o': BLAST tabular format.
            'sam': SAM format.
            'map': Simple map of query <tab> subject.
            unspecified: program will automatically infer from file content.
    input_ext : str, optional
        Input filename extension following sample Id.
            e.g., with input_ext = '_R1.fastq.gz', filename 'ID_R1.fastq.gz'
            is accepted and sample Id being 'ID', while filenames 'ID.log' and
            'readme.txt' are dropped.
    demux : bool, optional
        Whether demultiplex query Ids by pattern "sample_read".
            Query Id is split at first underscore. This works for QIIME/Qiita-
            style sample Ids (where underscore is not allowed). But one needs
            to be cautious otherwise.

    Returns
    -------
    dict of dict
        Per-rank profiles generated from classification

    See Also
    --------
    cli.classify

    Notes
    -----
    Explanations of parameters are provided as click decorators in `cli.py`.
    """
    # parse sample Ids
    if sample_ids:
        samples = read_ids(sample_ids)
        click.echo(f'Samples to include: {len(samples)}.')
    else:
        samples = None

    # match sample Ids and input alignment files
    demux, files = match_sample_file(input_path, input_ext, demux, samples)
    click.echo(f'Alignment files to read: {len(files)}.')

    # load gene coordinates
    if coords_fp:
        with readzip(coords_fp) as fh:
            coords = read_gene_coords(fh, sort=True)
        is_prefix = whether_prefix(coords)
    else:
        coords = None

    # load classification system
    args = (map_fps, names_fp, nodes_fp, newick_fp, ranktb_fp, lineage_fp)
    if any(args):
        click.echo(f'Reading classification system...', nl=False)
        tree, rankd, named, root = build_tree(*args)
        click.echo(' Done.')
        click.echo(f'Total classification units: {len(tree)}.')
    else:
        tree, rankd, named, root = None, None, None, None

    # parse target ranks
    ranks = ['none'] if rank_lst is None else rank_lst.split(',')
    data = {x: {} for x in ranks}

    # parse input maps and generate profile
    args = [tree, rankd, root]
    for fp in sorted(files):
        click.echo(f'Parsing alignment file {basename(fp)}...', nl=False)

        # read alignment file into query-subject(s) map
        with readzip(fp) as fh:

            # plain read to subject map
            if coords is None:
                readmap = read_align_file(fh, input_fmt)

            # match reads to genes using ordinal algorithm
            else:
                readmap = ordinal(fh, coords, input_fmt, prefix=is_prefix)

        click.echo(' Done.')
        click.echo(f'Query sequences: {len(readmap)}.')

        # merge duplicate subjects per query
        readmap = {k: set(v) for k, v in readmap.items()}

        # demultiplex into multiple samples
        if demux:
            readmap = demultiplex(readmap, samples)

        # sample Id from filename
        else:
            readmap = {files[fp]: readmap}

        # assign reads at each rank
        for sample, map_ in readmap.items():
            for rank in ranks:

                # call assignment workflow
                asgmt = {k: assign(v, rank, *args) for k, v in map_.items()}

                # convert floats into intergers
                counts = intize(count(asgmt))

                # delete "None" keys
                delnone(counts)

                # combine old and new counts
                try:
                    add_dict(data[rank][sample], counts)
                except KeyError:
                    data[rank][sample] = counts

    # get sample Ids from data
    if demux and not samples:
        samples = sorted(allkeys(data))

    # output results
    if output_path:

        # determine output filenames
        click.echo('Writing output profiles...', nl=False)
        if len(ranks) == 1:
            rank2fp = {ranks[0]: output_path}
        else:
            makedirs(output_path, exist_ok=True)
            rank2fp = {x: join(output_path, f'{x}.tsv') for x in ranks}

        # write output profile(s)
        for rank in ranks:
            with open(rank2fp[rank], 'w') as fh:
                write_profile(fh, data[rank], named, samples)
        click.echo(' Done.')

    click.echo('Task completed.')
    return data


def count(matches):
    """Count occurrences of taxa in a map.

    Parameters
    ----------
    matches : dict of str or dict
        Query-to-taxon(a) map.

    Returns
    -------
    dict
        Taxon-to-count map.
    """
    res = {}
    for taxa in matches.values():
        try:
            # unique match (scalar)
            res[taxa] = res.get(taxa, 0) + 1
        except TypeError:
            # multiple matches (dict of subject : count), to be normalized by
            # total match count
            k = 1 / sum(taxa.values())
            for taxon, n in taxa.items():
                res[taxon] = res.get(taxon, 0) + n * k
    return res


def majority(taxa, th=0.8):
    """Select taxon from list by majority rule.

    Parameters
    ----------
    taxa : list of str
        Input taxon list.
    th : float, optional
        Threshold of majority, range = (0.5, 1.0].

    Returns
    -------
    str or None
        Selected taxon.
    """
    for taxon, n in sorted(count_list(taxa).items(), key=lambda x: x[1],
                           reverse=True):
        return taxon if n >= len(taxa) * th else None


def assign(subs, rank=None, tree=None, rankd=None, root=None,
           above=False, multi=False, major=None, subok=False):
    """Assign a query sequence to a classification unit based on its subjects.

    Parameters
    ----------
    subs : set of str
        Subject(s) of a query sequence.
    rank : str, optional
        Ranks to assign to, or "free" for rank-free LCA assignment.
    tree : dict, optional
        Hierarchical classification system.
    rankd : dict, optional
        Rank dictionary.
    root : str, optional
        Root identifier.
    above : bool, optional
        Assignment above given rank is acceptable (for fixed ranks).
    major : float, optional
        Majority-rule assignment threshold (available only with a fixed rank
        and not above or multi).
    multi : bool, optional
        Count occurrence of each possible assignment instead of targeting one
        assignment (available only with a fixed rank and not above).
    subok : bool, optional
        Allow assignment to subject(s) itself instead of higher classification
        units.

    Returns
    -------
    str or dict
        Unique assignment or assignment-to-count map.
    """
    # no classification, just subject(s) itself
    if rank is None or rank == 'none' or tree is None:
        if len(subs) == 1:
            return max(subs)
        elif multi:
            return count_list(subs)
        else:
            return None

    # free rank classification: find LCA
    elif rank == 'free':
        if len(subs) == 1:
            if subok:
                return tree[max(subs)]
            else:
                return max(subs)
        else:
            lca = find_lca(subs, tree)
            return None if lca == root else lca

    # fixed rank classification
    else:
        taxa = [find_rank(x, rank, tree, rankd) for x in subs]
        if len(set(taxa)) == 1:
            return taxa[0]
        elif major:
            return majority(taxa, major)
        elif above:
            lca = find_lca(set(taxa), tree)
            return None if lca == root else lca
        elif multi:
            return count_list(taxa)
        else:
            return None


def write_profile(fh, data, named=None, samples=None):
    """Write profile to a plain tab-delimited file.
    """
    if samples is None:
        samples = sorted(data)
    print('#FeatureID\t{}'.format('\t'.join(samples)), file=fh)
    for key in sorted(allkeys(data)):
        # get feature name
        try:
            row = [named[key]]
        except (TypeError, KeyError):
            row = [key]
        # get feature count
        for sample in samples:
            try:
                row.append(str(data[sample][key]))
            except KeyError:
                row.append('0')
        print('\t'.join(row), file=fh)


def prep_table(profile, samples=None):
    """Convert a profile into data, index and columns, which can be further
    converted into a Pandas DataFrame or BIOM table.

    Parameters
    ----------
    profile : dict
        Input profile.

    Returns
    -------
    list of list, list, list
        Data (2D array of values).
        Index (observation Ids).
        Columns (sample Ids).
    """
    index = sorted(allkeys(profile))
    columns = samples or sorted(profile)
    data = []
    for key in index:
        row = []
        for sample in columns:
            try:
                row.append(profile[sample][key])
            except KeyError:
                row.append(0)
        data.append(row)
    return data, index, columns
