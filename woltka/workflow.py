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
import click

from .util import update_dict, allkeys, add_dict, intize, delnone
from .file import (
    readzip, path2stem, read_ids, id2file_from_dir, id2file_from_map,
    write_readmap, write_table)
from .align import Plain, parse_align_file
from .classify import (
    assign_none, assign_free, assign_rank, count, strip_index, demultiplex)
from .tree import (
    read_names, read_nodes, read_lineage, read_newick, read_ranktb, read_map,
    fill_root)
from .ordinal import Ordinal, read_gene_coords, whether_prefix
from .biom import profile_to_biom, write_biom


def workflow(input_fp:   str,
             output_fp:  str,
             # input
             input_fmt:    str = None,
             input_ext:    str = None,
             samples:      str = None,
             demux:       bool = None,
             lines:        int = 1000000,
             # hierarchies
             nodes_fp:     str = None,
             newick_fp:    str = None,
             lineage_fp:   str = None,
             ranktb_fp:    str = None,
             map_fps:     list = [],
             names_fp:     str = None,
             # assignment
             ranks:        str = None,
             above:       bool = False,
             major:       bool = None,
             ambig:       bool = True,
             subok:       bool = True,
             deidx:       bool = False,
             # gene matching
             coords_fp:    str = None,
             overlap:      int = 80,
             # output
             output_fmt:   str = None,
             name_as_id:  bool = False,
             add_rank:    bool = False,
             add_lineage: bool = False,
             outmap_dir:   str = None,
             outmap_zip:   str = 'gz') -> dict:
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

    # build classification system
    tree, rankdic, namedic, root = build_hierarchy(
        names_fp, nodes_fp, newick_fp, lineage_fp, ranktb_fp, map_fps)

    # build mapping module
    mapper = build_mapper(coords_fp, overlap)

    # target classification ranks
    ranks, rank2dir = prepare_ranks(ranks, outmap_dir)

    # classify query sequences
    data = classify(
        mapper, files, samples, input_fmt, demux, tree, rankdic, namedic, root,
        ranks, rank2dir, above, major and major / 100, ambig, subok, deidx,
        lines)

    # write output profiles
    write_profiles(
        data, output_fp, output_fmt, samples, tree, rankdic, namedic,
        name_as_id, add_rank, add_lineage)

    click.echo('Task completed.')
    return data


def classify(mapper:  object,
             files:     list or dict,
             samples:   list = None,
             input_fmt:  str = None,
             demux:     bool = None,
             tree:      dict = None,
             rankdic:   dict = None,
             namedic:   dict = None,
             root:       str = None,
             ranks:      str = None,
             rank2dir:  dict = None,
             above:     bool = False,
             major:      int = None,
             ambig:      str = True,
             subok:     bool = None,
             deidx:     bool = False,
             lines:      int = None) -> dict:
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
    input_fmt : str, optional
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

    Returns
    -------
    dict of dict
        Per-rank profiles generated from classification.
    """
    data = {x: {} for x in ranks}

    # assignment parameters
    kwargs = {'tree': tree, 'rankdic': rankdic, 'root':  root, 'above': above,
              'major': major, 'ambig': ambig, 'subok': subok,
              'namedic': namedic, 'rank2dir': rank2dir}

    # parse input maps and generate profile
    for fp in sorted(files):
        n = 0
        click.echo(f'Parsing alignment file {basename(fp)} ', nl=False)

        # read alignment file into query-subject(s) map
        with readzip(fp) as fh:

            # parse alignment file by chunk
            for rmap in parse_align_file(fh, mapper, input_fmt, lines):

                # show progress
                click.echo('.', nl=False)
                n += len(rmap)

                # reshape read map
                rmap = reshape_readmap(rmap, deidx, demux, samples, files, fp)

                # assign reads at each rank
                for sample, map_ in rmap.items():
                    for rank in ranks:

                        # call assignment workflow
                        assign_readmap(map_, data, rank, sample, **kwargs)

        click.echo(' Done.')
        click.echo(f'Number of query sequences: {n}.')

    click.echo('Task completed.')
    return data


def parse_samples(fp:      str,
                  ext:     str = None,
                  samples: str = None,
                  demux:  bool = None) -> (list, list or dict, bool):
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

    click.echo(f'Number of alignment files to read: {len(files)}.')
    click.echo(f'Demultiplexing: {"on" if demux else "off"}.')
    return samples, files, demux


def build_mapper(coords_fp: str = None,
                 overlap:   int = None) -> object:
    """Build mapping module (Plain or Ordinal).

    Parameters
    ----------
    coords_fp : str, optional
        Path to gene coordinates file.
    overlap : int, optional
        Read/gene overlapping percentage threshold.

    Returns
    -------
    object
        Mapping module.

    Notes
    -----
    Currently two mappers are supported: Plain() for regular alignments
    (i.e., simple query-to-subject maps), Ordinal() for alignments with
    coordinates which will be used to match queries (reads) and genes. The
    presence of a gene coordinates file (`coords_fp`) is an indicator for
    using the latter.
    """
    if coords_fp:
        click.echo(f'Reading gene coordinates...', nl=False)
        with readzip(coords_fp) as fh:
            coords = read_gene_coords(fh, sort=True)
        click.echo(' Done.')
        click.echo(f'Total number of host sequences: {len(coords)}.')
        return Ordinal(coords, whether_prefix(coords),
                       overlap and overlap / 100)
    else:
        return Plain()


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


def build_hierarchy(names_fp:   str = None,
                    nodes_fp:   str = None,
                    newick_fp:  str = None,
                    lineage_fp: str = None,
                    ranktb_fp:  str = None,
                    map_fps:   list = None) -> (dict, dict, dict, str):
    """Construct hierarchical classification system.

    Parameters
    ----------
    names_fp : str, optional
        Taxonomic names file.
    nodes_fp : str, optional
        Taxonomic nodes file.
    newick_fp : str, optional
        Newick tree file.
    lineage_fp : str, optional
        Lineage strings file.
    ranktb_fp : str, optional
        Rank table file.
    map_fps : list of str, optional
        Mapping files.

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
        names_fp, nodes_fp, newick_fp, lineage_fp, ranktb_fp, map_fps])
    if is_build:
        click.echo('Constructing classification system...')

    # taxonomy names
    if names_fp:
        click.echo(f'  Parsing taxonomy names file: {names_fp}...', nl=False)
        with readzip(names_fp) as f:
            namedic = read_names(f)
        click.echo(' Done.')

    # taxonomy nodes
    if nodes_fp:
        click.echo(f'  Parsing taxonomy nodes file: {nodes_fp}...', nl=False)
        with readzip(nodes_fp) as f:
            tree_, rankdic_ = read_nodes(f)
        update_dict(tree, tree_)
        update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # Newick-format tree
    if newick_fp:
        click.echo(f'  Parsing Newick tree file: {newick_fp}...', nl=False)
        with readzip(newick_fp) as f:
            update_dict(tree, read_newick(f))
        click.echo(' Done.')

    # lineage strings file
    if lineage_fp:
        click.echo(f'  Parsing lineage file: {lineage_fp}...', nl=False)
        with readzip(lineage_fp) as f:
            tree_, rankdic_ = read_lineage(f)
        update_dict(tree, tree_)
        update_dict(rankdic, rankdic_)
        click.echo(' Done.')

    # rank table file
    if ranktb_fp:
        click.echo(f'  Parsing rank table file: {ranktb_fp}...', nl=False)
        with readzip(ranktb_fp) as f:
            update_dict(tree, read_ranktb(f))
        click.echo(' Done.')

    # plain mapping files
    for fp in map_fps:
        click.echo(f'  Parsing simple map file: {fp}...', nl=False)
        rank = path2stem(fp)  # filename stem as rank
        with readzip(fp) as f:
            map_ = read_map(f)
        update_dict(tree, map_)
        # update_dict(rankdic, {k: rank for k in set(map_.values())})
        click.echo(' Done.')

    # fill root
    root = fill_root(tree)

    if is_build:
        click.echo('Classification system constructed.')
        click.echo(f'Total number of classification units: {len(tree)}.')
    return tree, rankdic, namedic, root


def reshape_readmap(rmap:    dict,
                    deidx:   bool = None,
                    demux:   bool = None,
                    samples: list = None,
                    files:   dict = None,
                    fp:       str = None) -> dict:
    """Reshape a read map.

    Parameters
    ----------
    rmap : dict
        Read map to manipulate.
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
        strip_index(rmap)

    # demultiplex into multiple samples
    if demux:
        return demultiplex(rmap, samples)

    # sample Id from filename
    else:
        return {files[fp]: rmap}


def assign_readmap(rmap:     dict,
                   data:     dict,
                   rank:     str,
                   sample:   str,
                   rank2dir: dict = None,
                   tree:     dict = None,
                   rankdic:  dict = None,
                   namedic:  dict = None,
                   root:      str = None,
                   above:    bool = False,
                   major:     int = None,
                   ambig:     str = True,
                   subok:    bool = None):
    """Assign query sequences in a query-to-subjects map to classification
    units based on their subjects.

    Parameters
    ----------
    rmap : dict
        Read map to assign.
    data : dict
        Master data structure.
    rank : str
        Target rank to assign to.
    sample : str
        Sample ID.
    rank2dir : dict, optional
        Directory of output maps per rank.
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
    """
    # determine assigner
    if rank is None or rank == 'none' or tree is None:
        assigner = assign_none
        args = (ambig,)
    elif rank == 'free':
        assigner = assign_free
        args = (tree, root, subok)
    else:
        assigner = assign_rank
        args = (rank, tree, rankdic, root, above, major, ambig)

    # call assigner
    asgmt = {k: assigner(v, *args) for k, v in rmap.items()}

    # write classification map
    if rank2dir is not None:
        with open(join(rank2dir[rank], f'{sample}.txt'), 'a') as f:
            write_readmap(f, asgmt, namedic)

    # count taxa
    counts = count(asgmt)

    # round floats
    intize(counts)

    # delete "None" keys
    delnone(counts)

    # combine old and new counts
    if sample in data[rank]:
        add_dict(data[rank][sample], counts)
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
        rank2fp = {x: join(fp, f'{x}.biom') for x in ranks}
        is_biom = is_biom is not False
    click.echo('Format of output feature table(s): {}.'.format(
        'BIOM' if is_biom else 'TSV'))

    # write output profile(s)
    click.echo('Writing output profiles...', nl=False)
    for rank, fp in rank2fp.items():
        if is_biom:
            write_biom(profile_to_biom(
                data[rank], samples, tree, rankdic, namedic), fp)
        else:
            with open(fp, 'w') as fh:
                write_table(fh, data[rank], samples,
                            tree if add_lineage else None,
                            rankdic if add_rank else None,
                            namedic, name_as_id)
    click.echo(' Done.')
