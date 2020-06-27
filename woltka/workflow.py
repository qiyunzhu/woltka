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
from collections import deque
from functools import partial, lru_cache
import click

from .util import update_dict, allkeys, sum_dict, intize
from .file import (
    openzip, path2stem, read_ids, id2file_from_dir, id2file_from_map, read_map,
    write_readmap, write_table)
from .align import plain_mapper
from .classify import (
    assign_none, assign_free, assign_rank, count, count_strata, strip_index,
    demultiplex)
from .tree import (
    read_names, read_nodes, read_lineage, read_newick, read_rank_table,
    fill_root)
from .ordinal import ordinal_mapper, read_gene_coords, whether_prefix
from .biom import profile_to_biom, write_biom


def workflow(input_fp:      str,
             output_fp:     str,
             # input
             input_fmt:     str = None,
             input_ext:     str = None,
             samples:       str = None,
             demux:        bool = None,
             lines:         int = 1000000,
             # hierarchies
             nodes_fp:      str = None,
             newick_fp:     str = None,
             lineage_fp:    str = None,
             rank_table_fp: str = None,
             map_fps:      list = [],
             map_as_rank:  bool = False,
             names_fps:    list = [],
             # assignment
             ranks:         str = None,
             above:        bool = False,
             major:        bool = None,
             ambig:        bool = True,
             subok:        bool = True,
             deidx:        bool = False,
             # gene matching
             coords_fp:     str = None,
             overlap:       int = 80,
             # stratification
             strata_dir:    str = None,
             # output
             output_fmt:    str = None,
             name_as_id:   bool = False,
             add_rank:     bool = False,
             add_lineage:  bool = False,
             outmap_dir:    str = None,
             outmap_zip:    str = 'gz') -> dict:
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
    # parse input samples
    samples, files, demux = parse_samples(
        input_fp, input_ext, samples, demux)

    # parse stratification files
    stratmap = parse_strata(strata_dir, samples)

    # build classification system
    tree, rankdic, namedic, root = build_hierarchy(
        names_fps, nodes_fp, newick_fp, lineage_fp, rank_table_fp, map_fps,
        map_as_rank)

    # build mapping module
    mapper = build_mapper(coords_fp, overlap)

    # target classification ranks
    ranks, rank2dir = prepare_ranks(ranks, outmap_dir)

    # classify query sequences
    data = classify(
        mapper, files, samples, input_fmt, demux, tree, rankdic, namedic, root,
        ranks, rank2dir, outmap_zip, above, major, ambig, subok, deidx, lines,
        stratmap)

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
             tree:      dict = None,
             rankdic:   dict = None,
             namedic:   dict = None,
             root:       str = None,
             ranks:      str = None,
             rank2dir:  dict = None,
             outzip:     str = None,
             above:     bool = False,
             major:      int = None,
             ambig:      str = True,
             subok:     bool = None,
             deidx:     bool = False,
             lines:      int = 1000000,
             stratmap:  dict = None) -> dict:
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
    deidx : bool, optional
        Strip "underscore index" suffixes from subject IDs.
    lines : int, optional
        Number of lines to read from alignment file per chunk.
    stratmap : dict, optional
        Map of sample ID to stratification file.

    Returns
    -------
    dict of dict
        Per-rank profiles generated from classification.
    """
    data = {x: {} for x in ranks}

    # assignment parameters
    kwargs = {'tree': tree, 'rankdic': rankdic, 'root':  root, 'above': above,
              'major': major and major / 100, 'ambig': ambig, 'subok': subok,
              'namedic': namedic, 'rank2dir': rank2dir, 'outzip': outzip if
              outzip != 'none' else None}

    # current sample Id
    csample = None

    # parse input maps and generate profile
    for fp in sorted(files):
        n = 0
        click.echo(f'Parsing alignment file {basename(fp)} ', nl=False)

        # read alignment file into query-subject(s) map
        with openzip(fp) as fh:

            # parse alignment file by chunk
            # for qryque, subque in parse_align_file(fh, mapper, fmt, lines):
            for qryque, subque in mapper(fh, fmt=fmt, n=lines):

                # show progress
                click.echo('.', nl=False)
                n += len(qryque)

                # reshape read map
                rmap = reshape_readmap(
                    qryque, subque, deidx, demux, samples, files, fp)

                # assign reads at each rank
                for sample, (qryque, subque) in rmap.items():

                    # read strata of current sample into cache
                    if stratmap and sample != csample:
                        with openzip(stratmap[sample]) as fh:
                            kwargs['strata'] = dict(read_map(fh))
                        csample = sample

                    for rank in ranks:

                        # call assignment workflow
                        assign_readmap(
                            qryque, subque, data, rank, sample, **kwargs)

        click.echo(' Done.')
        click.echo(f'  Number of query sequences: {n}.')

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
                 overlap:   int = None) -> callable:
    """Build mapper function (plain or ordinal).

    Parameters
    ----------
    coords_fp : str, optional
        Path to gene coordinates file.
    overlap : int, optional
        Read/gene overlapping percentage threshold.

    Returns
    -------
    callable
        Mapper function.

    Notes
    -----
    Currently two mappers are supported: "plain" for regular alignments
    (i.e., simple query-to-subject maps), "ordinal" for alignments with
    coordinates which will be used to match queries (reads) and genes. The
    presence of a gene coordinates file (`coords_fp`) is an indicator for
    using the latter.
    """
    if coords_fp:
        click.echo('Reading gene coordinates...', nl=False)
        with openzip(coords_fp) as fh:
            coords = read_gene_coords(fh, sort=True)
        click.echo(' Done.')
        click.echo(f'Total number of host sequences: {len(coords)}.')
        return partial(ordinal_mapper, coords=coords,
                       prefix=whether_prefix(coords),
                       th=overlap and overlap / 100)
    else:
        return plain_mapper


def prepare_ranks(ranks:      str = None,
                  outmap_dir: str = None) -> (list, dict or None):
    """Prepare classification ranks and directories of read-to-feature maps.

    Parameters
    ----------
    ranks : str, optional
        Target ranks (comma-separated).
    outmap_dir : str, optional
        Path to output read map directory.

    Returns
    -------
    list
        Ranks.
    dict or None
        Rank-to-directory map (or None if not necessary).
    """
    # determine ranks
    ranks = ['none'] if ranks is None else ranks.split(',')
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


def build_hierarchy(names_fps:    list = [],
                    nodes_fp:      str = None,
                    newick_fp:     str = None,
                    lineage_fp:    str = None,
                    rank_table_fp: str = None,
                    map_fps:      list = [],
                    map_as_rank:  bool = False) -> (dict, dict, dict, str):
    """Construct hierarchical classification system.

    Parameters
    ----------
    names_fps : list of str, optional
        Taxonomic names file(s).
    nodes_fp : str, optional
        Taxonomic nodes file.
    newick_fp : str, optional
        Newick tree file.
    lineage_fp : str, optional
        Lineage strings file.
    rank_table_fp : str, optional
        Rank table file.
    map_fps : list of str, optional
        Mapping file(s).
    map_as_rank : bool, optional
        Treat mapping filename stem as rank.

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
    is_build = any([
        names_fps, nodes_fp, newick_fp, lineage_fp, rank_table_fp, map_fps])
    if is_build:
        click.echo('Constructing classification system...')

    # taxonomy names
    for fp in names_fps:
        click.echo(f'  Parsing taxonomy names file: {fp}...', nl=False)
        with openzip(fp) as f:
            names = read_names(f)
        update_dict(namedic, names)
        click.echo(' Done.')

    # taxonomy nodes
    if nodes_fp:
        click.echo(f'  Parsing taxonomy nodes file: {nodes_fp}...', nl=False)
        with openzip(nodes_fp) as f:
            tree_, rankdic_ = read_nodes(f)
        update_dict(tree, tree_)
        update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # Newick-format tree
    if newick_fp:
        click.echo(f'  Parsing Newick tree file: {newick_fp}...', nl=False)
        with openzip(newick_fp) as f:
            update_dict(tree, read_newick(f))
        click.echo(' Done.')

    # lineage strings file
    if lineage_fp:
        click.echo(f'  Parsing lineage file: {lineage_fp}...', nl=False)
        with openzip(lineage_fp) as f:
            tree_, rankdic_ = read_lineage(f)
        update_dict(tree, tree_)
        update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # rank table file
    if rank_table_fp:
        click.echo(f'  Parsing rank table file: {rank_table_fp}...', nl=False)
        with openzip(rank_table_fp) as f:
            tree_, rankdic_ = read_rank_table(f)
            update_dict(tree, tree_)
            update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # plain mapping files
    for fp in map_fps:
        click.echo(f'  Parsing simple map file: {fp}...', nl=False)
        with openzip(fp) as f:
            map_ = dict(read_map(f, multi=False))
        update_dict(tree, map_)

        # filename stem as rank
        if map_as_rank:
            rank = path2stem(fp)
            update_dict(rankdic, {k: rank for k in set(map_.values())})
        click.echo(' Done.')

    # fill root
    root = fill_root(tree)

    if is_build:
        click.echo('Classification system constructed.')
        click.echo(f'Total number of classification units: {len(tree)}.')
    return tree, rankdic, namedic, root


def reshape_readmap(qryque: deque,
                    subque: deque,
                    deidx:   bool = None,
                    demux:   bool = None,
                    samples: list = None,
                    files:   dict = None,
                    fp:       str = None) -> dict:
    """Reshape a read map.

    Parameters
    ----------
    qryque : deque
        Query queue to manipulate.
    subque : deque
        Subject(s) queue to manipulate.
    deidx : bool, optional
        Strip suffixes from subject IDs.
    demux : bool, optional
        Whether perform demultiplexing.
    samples : list of str, optional
        Sample ID list to include.
    files : list or dict
        Map of filepaths to sample IDs.
    fp : str
        Path to current alignment file.

    Returns
    -------
    dict
        Reshaped read map.
    """
    # strip indices from subjects
    if deidx:
        subque = strip_index(subque)

    # demultiplex into multiple samples
    if demux:
        return demultiplex(qryque, subque, samples)

    # sample Id from filename
    else:
        return {files[fp] if files else None: (qryque, subque)}


def assign_readmap(qryque:  deque,
                   subque:  deque,
                   data:     dict,
                   rank:      str,
                   sample:    str,
                   rank2dir: dict = None,
                   outzip:    str = None,
                   tree:     dict = None,
                   rankdic:  dict = None,
                   namedic:  dict = None,
                   root:      str = None,
                   above:    bool = False,
                   major:   float = None,
                   ambig:     str = True,
                   subok:    bool = None,
                   strata:   dict = None):
    """Assign query sequences in a query-to-subjects map to classification
    units based on their subjects.

    Parameters
    ----------
    qryque : deque
        Query queue to assign.
    subque : deque
        Subject(s) queue for assignment.
    data : dict
        Master data structure.
    rank : str
        Target rank to assign to.
    sample : str
        Sample ID.
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
    strata : dict, optional
        Read-to-feature map for stratification.
    """
    # determine assigner
    if rank is None or rank == 'none' or tree is None:
        assigner = lru_cache(maxsize=128)(partial(
            assign_none, ambig=ambig))
    elif rank == 'free':
        assigner = lru_cache(maxsize=128)(partial(
            assign_free, tree=tree, root=root, subok=subok))
    else:
        assigner = lru_cache(maxsize=128)(partial(
            assign_rank, rank=rank, tree=tree, rankdic=rankdic, root=root,
            above=above, major=major, ambig=ambig))

    # call assigner
    asgmt = {}
    for query, subjects in zip(qryque, subque):
        # res = assigner(subjects, *args)
        res = assigner(frozenset(subjects))
        if res is not None:
            asgmt[query] = res

    # write classification map
    if rank2dir is not None:
        outfp = join(rank2dir[rank], f'{sample}.txt')
        with openzip(f'{outfp}.{outzip}' if outzip else outfp, 'at') as fh:
            write_readmap(fh, asgmt, namedic)

    # count taxa
    counts = count_strata(asgmt, strata) if strata else count(asgmt)

    # round floats
    intize(counts)

    # combine old and new counts
    if sample in data[rank]:
        sum_dict(data[rank][sample], counts)
    else:
        data[rank][sample] = counts


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
        Output BIOM instead of TSV format.
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
    click.echo('Writing output profiles...')
    for rank, fp in rank2fp.items():
        if is_biom:
            biom = profile_to_biom(
                data[rank], samples, tree if add_lineage else None, rankdic
                if add_rank else None, namedic, name_as_id)
            write_biom(biom, fp)
            m, n = biom.shape
        else:
            with open(fp, 'w') as fh:
                n, m = write_table(
                    fh, data[rank], samples, tree if add_lineage else None,
                    rankdic if add_rank else None, namedic, name_as_id)
        click.echo(f'  Rank: {rank}, samples: {n}, features: {m}.')

    click.echo('Profiles written.')
