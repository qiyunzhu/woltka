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
from .sample import read_ids, match_sample_file
from .align import read_align_file, strip_index, demultiplex
from .ordinal import read_gene_coords, whether_prefix, ordinal
from .tree import build_tree, find_rank, find_lca


def classify(input_path:  str,
             output_path: str,
             input_fmt:   str = None,
             input_ext:   str = None,
             sample_ids: list = None,
             demux:      bool = None,
             rank_lst:    str = None,
             above:      bool = False,
             major:       int = None,
             ambig:       str = True,
             subok:      bool = None,
             lca:        bool = True,
             deidx:      bool = False,
             coords_fp:   str = None,
             map_fps:    list = None,
             nodes_fp:    str = None,
             newick_fp:   str = None,
             ranktb_fp:   str = None,
             lineage_fp:  str = None,
             names_fp:    str = None) -> dict:
    """Main classification workflow.

    Parameters
    ----------
    input_path : str
        Path to input alignment file or directory of alignment files.
    output_path : str
        Path to output profile file or directory of profile files.
        Will be a file if there is one target rank, or a directory of multiple
        files starting with rank name if there are multiple target ranks.

    input_fmt : str, optional
        Format of input alignment file. Options:
        - 'b6o': BLAST tabular format.
        - 'sam': SAM format.
        - 'map': Simple map of query <tab> subject.
        If None, program will automatically infer from file content.
    input_ext : str, optional
        Input filename extension following sample ID.
        For example, with input_ext = '_R1.fastq.gz', filename 'ID_R1.fastq.gz'
        is accepted and sample ID being 'ID', while filenames 'ID.log' and
        'readme.txt' are dropped.
    sample_ids : iterable of str, optional
        List of sample IDs to be included.
    demux : bool, optional
        Whether demultiplex query IDs by pattern "sample_read".
        Query ID is split at first underscore. This works for QIIME/Qiita-style
        sample IDs (where underscore is not allowed). But one needs to be
        cautious otherwise.
        If None, program will determine demultiplex behavior based on the input
        path: turn on if it's a file; turn off if it's a directory.

    rank_lst: list of str, optional
        List of ranks at each of which sequences are to be classified. Can also
        be "none" to omit classification (simply report subject IDs) or "free"
        to perform free-rank classification (LCA of subjects regardless of rank
        will be reported).
    above : bool, optional
        Allow assigning to a classification unit higher than given rank.
    major : int, optional
        Perform majority-rule assignment based on this percentage threshold.
        Range: [51, 99].
    ambig : bool, optional
        Allow one sequence to be assigned to multiple classification units at
        the same rank. The profile will be normalized by the number of matches
        per query sequence.
    subok : bool, optional
        Allow directly reporting subject ID(s) if a sequence cannot be assigned
        to any higher classification unit.
    lca : bool, optional
        Perform lowest common ancestor (LCA) assignment for non-unique matches.
    deidx : bool, optional
        Strip "underscore index" suffixes from subject IDs.

    coords_fp : str, optional
        Path to table of gene coordinates on reference genomes, with which
        sequence-to-genome alignments will be translated into sequence-to-gene
        mappings. Essential for functional classification.

    map_fps : list of str, optional
        Paths to maps of lower classification units to higher classification
        units. Filename stem will be treated as the rank name.
        Examples include nucleotides to host genomes, sequence IDs to taxonomy
        IDs, gene families to pathways, etc.
    nodes_fp : str, optional
        Path to NCBI-style nodes.dmp or compatible formats, providing a map of
        taxon (1st column) to parent (2nd column) and rank (3rd column,
        optional).
    newick_fp : str, optional
        Path to Newick-format tree file which defines rank-free hierarchies of
        classification.
    ranktb_fp : str, optional
        Path to table of classification units at each rank (column).
    lineage_fp : str, optional
        Path to map of lineage strings (semicolon-delimited taxa from high to
        low). Can be Greengenes-style taxonomy where rank codes such as "k__"
        will be parsed, or code-free strings.
    names_fp : str, optional
        Path to NCBI-style names.dmp or plain taxon-to-name map.

    Returns
    -------
    dict of dict
        Per-rank profiles generated from classification.

    See Also
    --------
    .cli.classify
        Command-line parameters and help information.

    Notes
    -----
    This is the only function in the entire program which directly interface
    with the user via `click.echo`, except for raising errors.
    """
    # parse sample Ids
    samples = None
    if sample_ids:
        samples = read_ids(sample_ids)
        click.echo(f'Samples to include: {len(samples)}.')

    # match sample Ids and input alignment files
    demux, files = match_sample_file(input_path, input_ext, demux, samples)
    click.echo(f'Alignment files to read: {len(files)}.')

    # load gene coordinates
    coords = None
    if coords_fp:
        with readzip(coords_fp) as fh:
            coords = read_gene_coords(fh, sort=True)
        is_prefix = whether_prefix(coords)

    # load classification system
    tree = None
    args = (map_fps, names_fp, nodes_fp, newick_fp, ranktb_fp, lineage_fp)
    if any(args):
        click.echo(f'Reading classification system...', nl=False)
        tree, rankd, named, root = build_tree(*args)
        click.echo(' Done.')
        click.echo(f'Total classification units: {len(tree)}.')

    # parse target ranks
    ranks = ['none'] if rank_lst is None else rank_lst.split(',')
    data = {x: {} for x in ranks}

    # assignment parameters
    kwargs = {'tree':  tree,
              'rankd': None if tree is None else rankd,
              'root':  None if tree is None else root,
              'above': above,
              'major': None if major is None else major / 100,
              'ambig': ambig,
              'subok': subok}

    # parse input maps and generate profile
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

        # strip indices from subjects
        if deidx:
            strip_index(readmap)

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
                asgmt = {k: assign(v, rank, **kwargs) for k, v in map_.items()}

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
        kwargs = {'named': None if tree is None else named,
                  'samples': samples}
        for rank in ranks:
            with open(rank2fp[rank], 'w') as fh:
                write_profile(fh, data[rank], **kwargs)
        click.echo(' Done.')

    click.echo('Task completed.')
    return data


def assign(subs:    set,
           rank:    str = None,
           tree:   dict = None,
           rankd:  dict = None,
           root:    str = None,
           above:  bool = False,
           major: float = None,
           ambig:  bool = False,
           subok:  bool = False,
           deidx:  bool = False):
    """Assign a query sequence to a classification unit based on its subjects.

    Parameters
    ----------
    subs : set of str
        Subject(s) of a query sequence.
    rank : str, optional
        Rank to assign to, or "free" for rank-free LCA assignment.
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
        and not above or ambig).
    ambig : bool, optional
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
        elif ambig:
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
        elif ambig:
            return count_list(taxa)
        else:
            return None


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


def write_profile(fh, data, named=None, samples=None):
    """Write profile to a tab-delimited file.

    Parameters
    ----------
    fh : file handle
        Output file.
    data : dict
        Profile data.
    named : dict, optional
        Taxon name dictionary.
    samples : list, optional
        Ordered sample ID list.
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
