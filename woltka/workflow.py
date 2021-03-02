#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Main classification workflow.

Notes
-----
Only in this script can functions directly interface with the user by screen
output (via `click`) and file input/output, except for raising errors.
"""

from os import makedirs
from os.path import join, basename, isfile, isdir
from collections import deque, defaultdict
from itertools import compress
from functools import partial, lru_cache
import click

from .util import update_dict, allkeys, sum_dict, intize_dict
from .file import (
    openzip, readzip, path2stem, stem2rank, read_ids, id2file_from_dir,
    id2file_from_map, read_map_uniq, read_map_1st, write_readmap)
from .align import plain_mapper
from .classify import (
    assign_none, assign_free, assign_rank, count, count_strata)
from .tree import (
    read_names, read_nodes, read_lineage, read_newick, read_columns,
    fill_root)
from .ordinal import ordinal_mapper, read_gene_coords, whether_prefix
from .table import prep_table, write_table


def workflow(input_fp:     str,
             output_fp:    str,
             # input
             input_fmt:    str = None,
             input_ext:    str = None,
             samples:      str = None,
             demux:       bool = None,
             trimsub:      str = None,
             # hierarchies
             nodes_fps:   list = [],
             newick_fps:  list = [],
             lineage_fps: list = [],
             columns_fps: list = [],
             map_fps:     list = [],
             map_as_rank: bool = False,
             names_fps:   list = [],
             # assignment
             ranks:        str = None,
             uniq:        bool = False,
             major:       bool = None,
             above:       bool = False,
             subok:       bool = False,
             # gene matching
             coords_fp:    str = None,
             overlap:      int = 80,
             # stratification
             strata_dir:   str = None,
             # output
             output_fmt:   str = None,
             unassigned:  bool = False,
             name_as_id:  bool = False,
             add_rank:    bool = False,
             add_lineage: bool = False,
             outmap_dir:   str = None,
             outmap_zip:   str = 'gz',
             # performance
             chunk:        int = None,
             cache:        int = 1024,
             no_exe:      bool = False) -> dict:
    """Main classification workflow which accepts command-line arguments.

    Returns
    -------
    dict
        Resulting profile.

    Notes
    -----
    This function directly parses groups of command-line parameters to the
    following sub-workflows. Paramters are briefly explained in `cli.classify`
    and detailed in individual sub-workflows.

    See Also
    --------
    .cli.classify
        Command-line arguments and help information.
    """
    # available external compressors
    zippers = None if no_exe else {}

    # parse input samples
    samples, files, demux = parse_samples(
        input_fp, input_ext, samples, demux)

    # parse stratification files
    stratmap = parse_strata(strata_dir, samples)

    # build classification system
    tree, rankdic, namedic, root = build_hierarchy(
        names_fps, nodes_fps, newick_fps, lineage_fps, columns_fps, map_fps,
        map_as_rank, zippers)

    # build mapping module
    mapper, chunk = build_mapper(coords_fp, overlap, chunk, zippers)

    # target classification ranks
    ranks, rank2dir = prepare_ranks(ranks, outmap_dir, tree, rankdic)

    # classify query sequences
    data = classify(
        mapper, files, samples, input_fmt, demux, trimsub, tree, rankdic,
        namedic if name_as_id else None, root, ranks, rank2dir, outmap_zip,
        uniq, major, above, subok, unassigned, stratmap, chunk, cache, zippers)

    # write output profiles
    write_profiles(
        data, output_fp, output_fmt, samples, tree, rankdic, namedic,
        name_as_id, add_rank, add_lineage)

    click.echo('Task completed.')
    return data


def classify(mapper:  object,
             files:     list or dict,
             samples:   list = None,
             fmt:        str = None,
             demux:     bool = None,
             trimsub:    str = None,
             tree:      dict = None,
             rankdic:   dict = None,
             namedic:   dict = None,
             root:       str = None,
             ranks:      str = None,
             rank2dir:  dict = None,
             outzip:     str = None,
             uniq:      bool = False,
             major:      int = None,
             above:     bool = False,
             subok:     bool = False,
             unasgd:    bool = False,
             stratmap:  dict = None,
             chunk:      int = None,
             cache:      int = 1024,
             zippers:   dict = None) -> dict:
    """Core of the classification workflow.

    Parameters
    ----------
    mapper : object
        Mapping module (Plain or Ordinal).
    files : list or dict
        Paths to input alignment files, if multiplexed, or dictionary of file
        paths to sample IDs, if per-sample.
    samples : list of str, optional
        Sample ID list to include.
    fmt : str, optional
        Format of input alignment file. Options:
        - 'b6o': BLAST tabular format.
        - 'sam': SAM format.
        - 'map': Simple map of query <tab> subject.
        If None, program will automatically infer from file content.
    demux : bool, optional
        Whether perform demultiplexing.
    trimsub : str, optional
        Trim subject IDs at the last given delimiter.
    tree : dict, optional
        Taxonomic tree.
    rankdic : dict, optional
        Rank dictionary.
    namedic : dict, optional
        Taxon name dictionary.
    root : str, optional
        Root identifier.
    ranks: list of str, optional
        List of ranks at each of which sequences are to be classified. Can also
        be "none" to omit classification (simply report subject IDs) or "free"
        to perform free-rank classification (LCA of subjects regardless of rank
        will be reported).
    rank2dir : dict, otional
        Write classification map per rank to directory.
    outzip : str, optional
        Output read map compression method (gz, bz2, xz or None).
    uniq : bool, optional
        Assignment must be unique. Otherwise, report all possible assignments
        and normalize counts (for none- and fixed-rank assignments).
    major : int, optional
        In given-rank classification, perform majority-rule assignment based on
        this percentage threshold. Range: [51, 99].
    above : bool, optional
        Allow assigning to a classification unit higher than given rank.
    subok : bool, optional
        In free-rank classification, allow assigning sequences to their direct
        subjects instead of higher classification units, if applicable.
    unasgd : bool, optional
        Report unassigned sequences.
    stratmap : dict, optional
        Map of sample ID to stratification file.
    chunk : int, optional
        Number of lines per chunk to read from alignment file.
    cache : int, optional
        LRU cache size for classification results at each rank.
    zippers : dict, optional
        External compression programs.

    Returns
    -------
    dict of dict
        Per-rank profiles generated from classification.

    Notes
    -----
    Subject(s) of each query are converted into a frozenset. This is because
    frozenset is hashable, a property necessary for subsequent assignment
    result caching.
    """
    data = {x: {} for x in ranks}

    # assigners for each rank
    assigners = {}

    # assignment parameters
    kwargs = {'assigners': assigners, 'cache': cache, 'tree': tree, 'rankdic':
              rankdic, 'namedic': namedic, 'root':  root, 'uniq': uniq,
              'major': major and major / 100, 'above': above, 'subok': subok,
              'unasgd': unasgd, 'rank2dir': rank2dir, 'outzip': outzip if
              outzip != 'none' else None}

    # current sample Id
    csample = False

    # parse input alignment file(s) and generate profile(s)
    for fp in sorted(files):
        click.echo(f'Parsing alignment file {basename(fp)} ', nl=False)

        # read alignment file into query-to-subject(s) map
        with readzip(fp, zippers) as fh:

            # query and progress counters
            nqry, nstep = 0, -1

            # parse alignment file by chunk
            for qryque, subque in mapper(fh, fmt=fmt, n=chunk):
                nqry += len(qryque)

                # (optional) strip indices and freeze sets
                subque = deque(strip_suffix(subque, trimsub) if trimsub else
                               map(frozenset, subque))

                # (optional) demultiplex and generate per-sample maps
                rmaps = demultiplex(qryque, subque, samples) if demux else {
                    files[fp] if files else None: (qryque, subque)}

                # assign reads at each rank
                for sample, rmap in rmaps.items():

                    # (optional) read strata of current sample into cache
                    if stratmap and sample != csample:
                        kwargs['strata'] = read_strata(
                            stratmap[sample], zippers)
                        csample = sample

                    # call assignment workflow for each rank
                    for rank in ranks:
                        assign_readmap(*rmap, data, rank, sample, **kwargs)

                # show progress
                istep = nqry // 1000000 - nstep
                if istep:
                    click.echo('.' * istep, nl=False)
                    nstep += istep

        # round values and drop zeros
        round_profiles(data, uniq, major, above)

        click.echo(' Done.')
        click.echo(f'  Number of sequences classified: {nqry}.')

    click.echo('Classification completed.')
    return data


def parse_samples(fp:        str,
                  ext:       str = None,
                  samples:   str = None,
                  demux:    bool = None) -> (list, list or dict, bool):
    """Determine sample IDs, aligment files, and multiplex status.

    Parameters
    ----------
    fp : str
        Path to a file or a directory.
    ext : str, optional
        Filename extension.
    samples : str, optional
        Comma-separated list of sample IDs, or path to sample ID list file.
    demux : bool, optional
        Whether perform demultiplexing.

    Returns
    -------
    list
        Sample IDs to include.
    list or dict
        Filepaths if demultiplexing, or filepath to sample ID map if not.
    bool
        Whether perform demultiplexing.
    """
    # read sample Ids
    if samples:
        samples = read_ids(samples) if isfile(samples) else samples.split(',')
        click.echo(f'Number of samples to include: {len(samples)}.')

    errmsg = 'Provided sample IDs and actual files are inconsistent.'

    # path is a directory
    if isdir(fp):

        # turn off demultiplexing if not decided
        demux = demux or False

        # get a map of plausible sample Ids to files
        map_ = id2file_from_dir(fp, ext, not demux and samples)
        if len(map_) == 0:
            raise ValueError('No valid file found in directory.')
        if demux:
            files = sorted([join(fp, x) for x in map_.values()])

        # validate with given sample Ids
        else:
            if not samples:
                samples = sorted(map_.keys())
            elif len(map_) < len(samples):
                raise ValueError(errmsg)
            files = {join(fp, map_[x]): x for x in samples}

        click.echo(f'Input directory: {fp}.')
        click.echo(f'Number of alignment files to read: {len(files)}.')

    # path is a file
    elif isfile(fp):

        # check if file is an Id-to-file map
        map_ = id2file_from_map(fp)
        if map_:

            # turn off demultiplexing if not decided
            demux = demux or False

            # validate with given sample Ids
            if samples:
                map_ = dict(map_)
                try:
                    files = {map_[x]: x for x in samples}
                except KeyError:
                    raise ValueError(errmsg)

            # sample order and files provided
            else:
                samples = [x[0] for x in map_]
                files = {x[1]: x[0] for x in map_}

            click.echo(f'Number of alignment files to read: {len(files)}.')

        # treat file as a single alignment file
        else:

            # turn on demultiplexing if not decided
            demux = demux is not False
            if demux:
                files = [fp]

            # validate with given sample Ids
            else:
                sample = path2stem(fp, ext)
                if samples and samples != [sample]:
                    raise ValueError(errmsg)
                files = {fp: sample}
                samples = [sample]

            click.echo(f'Input alignment file: {fp}.')

    else:
        raise ValueError(f'"{fp}" is not a valid file or directory.')

    click.echo(f'Demultiplexing: {"on" if demux else "off"}.')

    return samples, files, demux


def parse_strata(fp:       str = None,
                 samples: list = None) -> dict:
    """Get a map of sample Ids to mapping files for stratification.

    Parameters
    ----------
    fp : str, optional
        Directory of read-to-feature maps for stratification.

    Returns
    -------
    dict or None
        Sample ID to stratification file map, or None if not applicable.

    Raises
    ------
    ValueError
        Sample IDs are given but stratification files are missing for one or
        more samples.
    """
    if not fp:
        return
    click.echo(f'Stratification file directory: {fp}.')
    map_ = id2file_from_dir(fp, ids=samples)
    if len(samples or []) > len(map_):
        raise ValueError(
            'Cannot locate stratification files for one or more samples.')
    return {k: join(fp, v) for k, v in map_.items()}


def build_mapper(coords_fp: str = None,
                 overlap:   int = None,
                 chunk:     int = None,
                 zippers:  dict = None) -> (callable, int):
    """Build mapper function (plain or ordinal).

    Parameters
    ----------
    coords_fp : str, optional
        Path to gene coordinates file.
    overlap : int, optional
        Read/gene overlapping percentage threshold.
    chunk : int, optional
        Number of lines per chunk to read from alignment file.
    zippers : dict, optional
        External compression programs.

    Returns
    -------
    callable
        Mapper function.
    int
        Number of lines per chunk.

    Notes
    -----
    Currently two mappers are supported: "plain" for regular alignments
    (i.e., simple query-to-subject maps), "ordinal" for alignments with
    coordinates which will be used to match queries (reads) and genes. The
    presence of a gene coordinates file (`coords_fp`) is an indicator for
    using the latter.

    The chunk size, if not specified, is determined empirically: 1,000 for
    plain mapper and 1,000,000 for ordinal mapper.
    """
    if coords_fp:
        click.echo('Reading gene coordinates...', nl=False)
        with readzip(coords_fp, zippers) as fh:
            coords = read_gene_coords(fh, sort=True)
        click.echo(' Done.')
        click.echo(f'  Total number of host sequences: {len(coords)}.')
        chunk = chunk or 1000000
        return partial(ordinal_mapper, coords=coords,
                       prefix=whether_prefix(coords),
                       th=overlap and overlap / 100), chunk
    else:
        chunk = chunk or 1000
        return plain_mapper, chunk


def prepare_ranks(ranks:      str = None,
                  outmap_dir: str = None,
                  tree:      dict = None,
                  rankdic:   dict = None) -> (list, dict or None):
    """Prepare classification ranks and directories of read-to-feature maps.

    Parameters
    ----------
    ranks : str, optional
        Target ranks (comma-separated).
    outmap_dir : str, optional
        Path to output read map directory.
    tree : dict, optional
        Taxonomic tree.
    rankdic : dict, optional
        Rank dictionary.

    Returns
    -------
    list
        Ranks.
    dict or None
        Rank-to-directory map (or None if not necessary).
    """
    # ranks are given
    if ranks:
        ranks = ranks.split(',')

        # check presence of ranks in classification system
        if rankdic is not None:
            missing = set(ranks) - {'free', 'none'} - set(rankdic.values())
            if missing:
                raise ValueError(f'Ranks {", ".join(sorted(missing))} are not'
                                 ' found in classification system.')

    # if classification system is provided, do free-rank classification;
    # otherwise do no-rank assignment.
    else:
        ranks = ['free' if tree else 'none']

    click.echo('Classification will operate on these ranks: {}.'.format(
        ', '.join(ranks)))

    # check output directory
    if not outmap_dir:
        return ranks, None
    makedirs(outmap_dir, exist_ok=True)
    click.echo(f'Read-to-feature maps will be saved to: {outmap_dir}.')

    # determine output read map directory per rank
    if len(ranks) == 1:
        rank2dir = {ranks[0]: outmap_dir}
    else:
        rank2dir = {}
        for rank in ranks:
            dir_ = join(outmap_dir, rank)
            makedirs(dir_, exist_ok=True)
            rank2dir[rank] = dir_

    return ranks, rank2dir


def build_hierarchy(names_fps:   list = [],
                    nodes_fps:   list = [],
                    newick_fps:  list = [],
                    lineage_fps: list = [],
                    columns_fps: list = [],
                    map_fps:     list = [],
                    map_as_rank: bool = False,
                    zippers:     dict = None) -> (dict, dict, dict, str):
    """Construct hierarchical classification system.

    Parameters
    ----------
    names_fps : list of str, optional
        Taxonomic names file(s).
    nodes_fps : list of str, optional
        Taxonomic nodes file(s).
    newick_fps : list of str, optional
        Newick tree file.
    lineage_fps : list of str, optional
        Lineage strings file.
    columns_fps : list of str, optional
        Rank-per-column file.
    map_fps : list of str, optional
        Mapping file(s).
    map_as_rank : bool, optional
        Treat mapping filename stem as rank.
    zippers : dict, optional
        External compression programs.

    Returns
    -------
    dict
        Taxonomic tree.
    dict
        Rank dictionary.
    dict
        Name dictionary.
    str
        Root identifier.
    """
    tree, rankdic, namedic = {}, {}, {}

    # check if at least one filepath is specified
    is_build = any([
        names_fps, nodes_fps, newick_fps, lineage_fps, columns_fps, map_fps])
    if is_build:
        click.echo('Constructing classification system...')

    # taxonomy names
    for fp in names_fps:
        click.echo(f'  Parsing taxon names file: {basename(fp)}...', nl=False)
        with readzip(fp, zippers) as f:
            names = read_names(f)
        update_dict(namedic, names)
        click.echo(' Done.')

    # taxonomy nodes
    for fp in nodes_fps:
        click.echo(f'  Parsing taxon nodes file: {basename(fp)}...', nl=False)
        with readzip(fp, zippers) as f:
            tree_, rankdic_ = read_nodes(f)
        update_dict(tree, tree_)
        update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # Newick-format tree
    for fp in newick_fps:
        click.echo(f'  Parsing Newick tree file: {basename(fp)}...', nl=False)
        with readzip(fp, zippers) as f:
            update_dict(tree, read_newick(f))
        click.echo(' Done.')

    # lineage file
    for fp in lineage_fps:
        click.echo(f'  Parsing lineage file: {basename(fp)}...', nl=False)
        with readzip(fp, zippers) as f:
            tree_, rankdic_ = read_lineage(f)
        update_dict(tree, tree_)
        update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # columns file
    for fp in columns_fps:
        click.echo(f'  Parsing columns file: {basename(fp)}...', nl=False)
        with readzip(fp, zippers) as f:
            tree_, rankdic_ = read_columns(f)
            update_dict(tree, tree_)
            update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # plain mapping files
    for fp in map_fps:
        click.echo(f'  Parsing simple map file: {basename(fp)}...', nl=False)
        with readzip(fp, zippers) as f:
            map_ = dict(read_map_1st(f))
        update_dict(tree, map_)

        # filename stem as rank
        if map_as_rank:
            rank = stem2rank(path2stem(fp))
            update_dict(rankdic, {k: rank for k in set(map_.values())})
        click.echo(' Done.')

    # fill root
    root = fill_root(tree)

    if is_build:
        click.echo('Classification system constructed.')
        click.echo(f'  Total number of classification units: {len(tree)}.')

    return tree, rankdic, namedic, root


def strip_suffix(subque: list,
                 sep:     str) -> object:
    """Strip suffixes from subject IDs at the last separator.

    Parameters
    ----------
    subque : iterable
        Subject(s) queue to manipulate.
    sep : str, optional
        Separator between subject ID and suffix.

    Returns
    -------
    generator of frozenset
        Processed subject(s) queue.

    Notes
    -----
    This function will find the last occurrence of separator in a subject Id,
    and trim from it to the right end. If not found, the whole subject Id will
    be retained.
    """
    return map(frozenset, map(partial(
        map, lambda x: x.rsplit(sep, 1)[0]), subque))


def demultiplex(qryque:  list,
                subque:  list,
                samples: list = None,
                sep:      str = '_') -> dict:
    """Demultiplex a query-to-subject(s) map.

    Parameters
    ----------
    qryque : iterable
        Query queue to demultiplex.
    subque : iterable
        Corresponding subject(s) queue.
    samples : iterable of str, optional
        Sample IDs to keep.
    sep : str, optional
        Separator between sample ID and read ID.

    Returns
    -------
    dict of (deque, deque)
        Per-sample read-to-subject(s) maps.

    Notes
    -----
    In a multiplexed alignment file, query IDs are composed of sample ID and
    read ID, separated by a character (default: "_"). This function separates
    them at the first occurrence of the separator from left. If the separator
    is not found, the entire query ID will be retained as read ID and sample
    ID will be `None`.
    """
    if samples:
        samset = set(samples)

    # per-sample read and subject(s) queues
    qryques, subques = defaultdict(deque), defaultdict(deque)

    # current sample Id (it can be None so start with False)
    csample = False

    # list append method references
    qry_add, sub_add = None, None

    for query, subjects in zip(qryque, subque):

        # split query Id by first separator
        left, _, right = query.partition(sep)

        # if separator is present, take left and right
        # if there is no separator, take None and left
        sample, read = right and left, right or left

        # append read Id and subject(s) to queues
        if sample == csample:
            qry_add(read)
            sub_add(subjects)

        # (re-)assign method references to current sample
        elif not samples or sample in samset:
            csample = sample
            qry_add = qryques[sample].append
            sub_add = subques[sample].append

            qry_add(read)
            sub_add(subjects)

    return {x: (qryques[x], subques[x]) for x in qryques}


def read_strata(strata_fp: str,
                zippers:  dict = None) -> dict:
    """Read a stratification file for a sample.

    Parameters
    ----------
    strata_fp : str,
        Path to stratification file.
    zippers : dict, optional
        External compression programs.

    Returns
    -------
    dict
        Stratification information.

    Raises
    ------
    ValueError
        No stratification information is found in file.
    """
    with readzip(strata_fp, zippers) as fhs:
        strata = dict(read_map_uniq(fhs))
    if not strata:
        raise ValueError('No stratification information is found in file: '
                         f'{basename(strata_fp)}.')
    return strata


def assign_readmap(qryque:    list,
                   subque:    list,
                   data:      dict,
                   rank:       str,
                   sample:     str,
                   assigners: dict,
                   cache:      int = 1024,
                   rank2dir:  dict = None,
                   outzip:     str = None,
                   tree:      dict = None,
                   rankdic:   dict = None,
                   namedic:   dict = None,
                   root:       str = None,
                   uniq:      bool = False,
                   major:    float = None,
                   above:     bool = False,
                   subok:     bool = False,
                   unasgd:    bool = False,
                   strata:    dict = None):
    """Assign query sequences in a query-to-subjects map to classification
    units based on their subjects.

    Parameters
    ----------
    qryque : iterable of str
        Query queue to assign.
    subque : iterable of frozenset
        Subject(s) queue for assignment.
    data : dict
        Master data structure.
    rank : str
        Target rank to assign to.
    sample : str
        Sample ID.
    assigners : dict of callable
        Per-rank assigners.
    cache : int, optional
        LRU cache size for classification results at each rank.
    rank2dir : dict, optional
        Directory of output maps per rank.
    outzip : str, optional
        Output read map compression method (gz, bz2, xz or None).
    tree : dict, optional
        Hierarchical classification system.
    rankdic : dict, optional
        Rank dictionary.
    namedic : dict
        Taxon name directory.
    root : str, optional
        Root identifier.
    uniq : bool, optional
        Assignment must be unique. Otherwise, report all possible assignments
        and normalize counts (for none- and fixed-rank assignments).
    major : float, optional
        In given-rank classification, perform majority-rule assignment based on
        this fraction threshold.
    above : bool, optional
        In given-rank classification, assignment above the specified rank is
        acceptable.
    subok : bool, optional
        In free-rank classification, allow assigning sequences to their direct
        subjects instead of higher classification units, if applicable.
    unasgd : bool, optional
        Report unassigned sequences.
    strata : dict, optional
        Read-to-feature map for stratification.
    """
    # determine assigner and initiate (if not already)
    if rank is None or rank == 'none' or tree is None:
        if 'none' not in assigners:
            assigners['none'] = lru_cache(maxsize=cache)(partial(
                assign_none, uniq=uniq))
        assigner = assigners['none']
    elif rank == 'free':
        if 'free' not in assigners:
            assigners['free'] = lru_cache(maxsize=cache)(partial(
                assign_free, tree=tree, root=root, subok=subok))
        assigner = assigners['free']
    else:
        if rank not in assigners:
            assigners[rank] = lru_cache(maxsize=cache)(partial(
                assign_rank, rank=rank, tree=tree, rankdic=rankdic, root=root,
                major=major, above=above, uniq=uniq))
        assigner = assigners[rank]

    # call assigner on suject(s) per query
    resque = map(assigner, subque)

    # report or drop unassigned
    if unasgd:
        resque = [x or 'Unassigned' for x in resque]
    else:
        resque = list(resque)
        keep = list(map(None.__ne__, resque))
        qryque, resque = list(compress(qryque, keep)), list(compress(
            resque, keep))

    # write classification map
    if rank2dir is not None:
        outfp = join(rank2dir[rank], f'{sample}.txt')
        with openzip(f'{outfp}.{outzip}' if outzip else outfp, 'at') as fh:
            write_readmap(fh, qryque, resque, namedic)

    # count taxa
    counts = count_strata(qryque, resque, strata) if strata else count(resque)

    # combine old and new counts
    if sample in data[rank]:
        sum_dict(data[rank][sample], counts)
    else:
        data[rank][sample] = counts


def round_profiles(data:   list,
                   uniq:   bool = False,
                   major: float = None,
                   above:  bool = False):
    """Round counts in profiles into integers, and drop zero counts.

    Parameters
    ----------
    data : dict
        Profile data.
    uniq : bool, optional
        Assignment is unique.
    major : int, optional
        Majority-rule assignment threshold.
    above : bool, optional
        Assignment above given rank.
    """
    for rank in data:

        # free-rank classification is always unique
        if rank == 'free':
            continue

        # non-unique no-rank classification needs to be rounded
        elif rank == 'none':
            if uniq:
                continue
            for datum in data[rank].values():
                intize_dict(datum)

        # given rank classification needs rounding when all three criteria
        # are not met
        elif not any((uniq, major, above)):
            for datum in data[rank].values():
                intize_dict(datum)


def write_profiles(data:        dict,
                   fp:           str,
                   is_biom:     bool = None,
                   samples:     list = None,
                   tree:        dict = None,
                   rankdic:     dict = None,
                   namedic:     dict = None,
                   name_as_id:  bool = False,
                   add_rank:    bool = False,
                   add_lineage: bool = False):
    """Write profile to an output table file.

    Parameters
    ----------
    data : dict
        Profile data.
    fp : str
        Path to output file or directory.
    is_biom : bool, optional
        Output BIOM or TSV format.
    samples : list, optional
        Ordered sample ID list.
    tree : dict, optional
        Taxonomic tree.
    rankdic : dict, optional
        Rank dictionary.
    namedic : dict, optional
        Taxon name dictionary.
    name_as_id : bool, optional
        Replace feature IDs with names.
    add_lineage: bool optional
        Append feature ranks to table.
    add_lineage: bool optional
        Append lineage strings to table.

    Notes
    -----
    The boolean parameter `is_biom` can be one of the three values:
    - `None` (default): To be auto-determined based on the user-supplied output
    filename extension. If ".biom", do BIOM, otherwise do TSV. If output path
    is a directory, do BIOM by default.
    - `True` (command-line flag `--to-biom`): BIOM format.
    - `False` (command-line flag `--to-tsv`): TSV format.
    """
    if not fp:
        return

    # determine sample order
    if not samples:
        samples = sorted(allkeys(data))

    # determine ranks
    ranks = sorted(data.keys())

    # determine output filename and format
    if len(ranks) == 1:
        # single output file
        rank2fp = {ranks[0]: fp}
        if is_biom is None:
            is_biom = fp.endswith('.biom')
    else:
        # multiple output files
        makedirs(fp, exist_ok=True)
        is_biom = is_biom is not False
        ext = 'biom' if is_biom else 'tsv'
        rank2fp = {x: join(fp, f'{x}.{ext}') for x in ranks}
    click.echo('Format of output feature table(s): {}.'.format(
        'BIOM' if is_biom else 'TSV'))

    # override Id-to-name conversion
    if namedic is None:
        name_as_id = False

    # write output profile(s)
    fmt = 'BIOM' if is_biom else 'TSV'
    click.echo(f'Writing output profiles in {fmt} format...')
    for rank, fp in rank2fp.items():
        table = prep_table(data[rank], samples, tree if add_lineage else None,
                           rankdic if add_rank else None, namedic, name_as_id)
        write_table(table, fp, is_biom)
        n, m = len(table[2]), len(table[1])
        click.echo(f'  Rank: {rank}, samples: {n}, features: {m}.')

    click.echo('Profiles written.')
