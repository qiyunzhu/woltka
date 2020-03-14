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

from .util import (
    readzip, path2stem, update_dict, allkeys, read_ids, add_dict, intize,
    delnone, id2file_map, write_map)
from .align import Plain, parse_align_file
from .classify import (
    assign_none, assign_free, assign_rank, count, strip_index, demultiplex)
from .tree import (
    read_names, read_nodes, read_lineage, read_newick, read_ranktb, read_map,
    fill_root)
from .ordinal import Ordinal, read_gene_coords, whether_prefix
from .biom import profile_to_biom, write_biom


def workflow(input_path, output_path, outmap_dir=None,
             input_fmt=None, input_ext=None, sample_ids=None, demux=None,
             ranks=None, above=False, major=None, ambig=True, subok=True,
             deidx=False, coords_fp=None, overlap=80, names_fp=None,
             nodes_fp=None, newick_fp=None, lineage_fp=None, ranktb_fp=None,
             map_fps=[], map_rank=False, lines=1000000):
    """Main classification workflow which accepts command-line arguments.

    See Also
    --------
    .cli.classify
        Command-line arguments and help information.
    """
    # parse input samples
    samples, files, demux = parse_samples(
        input_path, input_ext, sample_ids, demux)

    # build classification system
    tree, rankd, named, root = build_hierarchy(
        names_fp, nodes_fp, newick_fp, lineage_fp, ranktb_fp, map_fps,
        map_rank)

    # build processor for alignment
    proc = build_align_proc(coords_fp, overlap)

    # target classification ranks
    ranks, rank2dir = prep_ranks(ranks, outmap_dir)

    # classify query sequences
    data = classify(
        proc, files, samples, input_fmt, demux, tree, rankd, named, root,
        ranks, rank2dir, above, major and major / 100, ambig, subok, deidx,
        lines)

    # write output profiles
    write_profiles(output_path, data, named, samples)
    click.echo('Task completed.')
    return data


def classify(proc:    object,
             files:     list or dict,
             samples:   list = None,
             input_fmt:  str = None,
             demux:     bool = None,
             tree:      dict = None,
             rankd:     dict = None,
             named:     dict = None,
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
    proc : object
        Alignment processing module.
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
    rankd : dict, optional
        Rank dictionary.
    named : dict, optional
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
    kwargs = {'tree':   tree, 'rankd': rankd,  'root':  root, 'above': above,
              'major': major, 'ambig': ambig, 'subok': subok, 'named': named,
              'rank2dir': rank2dir}

    # parse input maps and generate profile
    for fp in sorted(files):
        n = 0
        click.echo(f'Parsing alignment file {basename(fp)} ', nl=False)

        # read alignment file into query-subject(s) map
        with readzip(fp) as fh:

            # parse alignment file by chunk
            for rmap in parse_align_file(fh, proc, input_fmt, lines):
                click.echo('.', nl=False)
                n += len(rmap)

                # reshape read map
                rmap = reshape_rmap(rmap, deidx, demux, samples, files, fp)

                # assign reads at each rank
                for sample, map_ in rmap.items():
                    for rank in ranks:

                        # call assignment workflow
                        assign_rmap(map_, data, rank, sample, **kwargs)

        click.echo(' Done.')
        click.echo(f'Query sequences: {n}.')

    click.echo('Task completed.')
    return data


def parse_samples(path_:  str,
                  ext:    str = None,
                  ids:   list = None,
                  demux: bool = None) -> (list, list or dict, bool):
    """Determine sample IDs, aligment files, and multiplex status.

    Parameters
    ----------
    path_ : str
        Path to a file or a directory.
    ext : str, optional
        Filename extension.
    ids : iterable of str, optional
        List of sample IDs.
    demux : bool, optional
        Whether perform demultiplexing.

    Returns
    -------
    list, list or dict, bool
        Sample IDs to include.
        Filepaths if demultiplexing, or filepath to sample ID map if not.
        Whether perform demultiplexing.
    """
    # read sample Ids
    if ids:
        samples = read_ids(ids)
        click.echo(f'Samples to include: {len(samples)}.')
    else:
        samples = None

    errmsg = 'Provided sample IDs and actual files are inconsistent.'

    # path is an alignment file
    if isfile(path_):

        # turn on demultiplexing if not decided
        demux = demux is not False
        if demux:
            files = [path_]

        # validate with given sample Ids
        else:
            sample = path2stem(path_, ext)
            if samples and samples != [sample]:
                raise ValueError(errmsg)
            files = {path_: sample}
            samples = [sample]

    # path is a directory of alignment files
    elif isdir(path_):

        # turn off demultiplexing if not decided
        demux = demux or False

        # get a map of plausible sample Ids to files
        map_ = id2file_map(path_, ext, not demux and samples)
        if len(map_) == 0:
            raise ValueError('No valid file found in directory.')
        if demux:
            files = sorted([join(path_, x) for x in map_.values()])

        # validate with given sample Ids
        else:
            if not samples:
                samples = sorted(map_.keys())
            elif len(map_) < len(samples):
                raise ValueError(errmsg)
            files = {join(path_, map_[x]): x for x in samples}

    else:
        raise ValueError(f'"{path_}" is not a valid file or directory.')

    click.echo(f'Alignment files to read: {len(files)}.')
    click.echo(f'Demultiplexing: {"on" if demux else "off"}.')
    return samples, files, demux


def build_align_proc(coords_fp: str = None,
                     overlap:   int = None) -> object:
    """Build alignment processor.

    Parameters
    ----------
    coords_fp : str, optional
        Path to gene coordinates file.
    overlap : int, optional
        Read/gene overlapping percentage threshold

    Returns
    -------
    object
        Alignment processor.

    Notes
    -----
    Currently two processors are supported: Plain() for regular alignments
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
        click.echo(f'Total host sequences: {len(coords)}.')
        return Ordinal(coords, whether_prefix(coords),
                       overlap and overlap / 100)
    else:
        return Plain()


def build_hierarchy(names_fp:   str = None,
                    nodes_fp:   str = None,
                    newick_fp:  str = None,
                    lineage_fp: str = None,
                    ranktb_fp:  str = None,
                    map_fps:   list = None,
                    map_rank:  bool = False) -> (dict, dict, dict, str):
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
    map_rank : bool, optional
        Mapping filename is rank.

    Returns
    -------
    dict, dict, dict, str
        Taxonomic tree.
        Rank dictionary.
        Name dictionary.
        Root identifier.
    """
    tree, rankd, named = {}, {}, {}
    is_build = any([
        names_fp, nodes_fp, newick_fp, lineage_fp, ranktb_fp, map_fps])
    if is_build:
        click.echo(f'Reading classification system...', nl=False)

    # taxonomy names
    if names_fp:
        with readzip(names_fp) as f:
            named = read_names(f)

    # taxonomy nodes
    if nodes_fp:
        with readzip(nodes_fp) as f:
            tree_, rankd_ = read_nodes(f)
        update_dict(tree, tree_)
        update_dict(rankd, rankd_)

    # Newick-format tree
    if newick_fp:
        with readzip(newick_fp) as f:
            update_dict(tree, read_newick(f))

    # lineage strings file
    if lineage_fp:
        with readzip(lineage_fp) as f:
            tree_, rankd_ = read_lineage(f)
        update_dict(tree, tree_)
        update_dict(rankd, rankd_)

    # rank table file
    if ranktb_fp:
        with readzip(ranktb_fp) as f:
            update_dict(tree, read_ranktb(f))

    # plain mapping files
    for fp in map_fps:
        rank = path2stem(fp)  # filename stem as rank
        with readzip(fp) as f:
            map_ = read_map(f)
        update_dict(tree, map_)
        if map_rank:
            update_dict(rankd, {k: rank for k in set(map_.values())})

    # fill root
    root = fill_root(tree)

    if is_build:
        click.echo(' Done.')
        click.echo(f'Total classification units: {len(tree)}.')
    return tree, rankd, named, root


def reshape_rmap(rmap:    dict,
                 deidx:   bool = None,
                 demux:   bool = None,
                 samples: list = None,
                 files:   dict = None,
                 fp:      str = None) -> dict:
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


def assign_rmap(rmap:     dict,
                data:     dict,
                rank:     str,
                sample:   str,
                rank2dir: dict = None,
                tree:     dict = None,
                rankd:    dict = None,
                named:    dict = None,
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
    rankd : dict, optional
        Rank dictionary.
    named : dict
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
        args = (rank, tree, rankd, root, above, major, ambig)

    # call assigner
    asgmt = {k: assigner(v, *args) for k, v in rmap.items()}

    # write classification map
    if rank2dir is not None:
        with open(join(rank2dir[rank], f'{sample}.tsv'), 'a') as f:
            write_map(f, asgmt, named)

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


def write_profiles(path_:    str,
                   data:    dict,
                   samples: list = None,
                   tree:    dict = None,
                   rankd:   dict = None,
                   named:   dict = None):
    """Write profile to an output file.

    Parameters
    ----------
    data : dict
        Profile data.
    path_ : str
        Path to output file or directory.
    samples : list, optional
        Ordered sample ID list.
    tree : dict, optional
        Taxonomic tree.
    rankd : dict, optional
        Rank dictionary.
    named : dict, optional
        Taxon name dictionary.
    """
    if not path_:
        return
    if not samples:
        samples = sorted(allkeys(data))
    ranks = sorted(data.keys())

    # determine output filenames
    click.echo('Writing output profiles...', nl=False)
    if len(ranks) == 1:
        rank2fp = {ranks[0]: path_}
    else:
        makedirs(path_, exist_ok=True)
        rank2fp = {x: join(path_, f'{x}.biom') for x in ranks}

    # write output profile(s)
    for rank, fp in rank2fp.items():
        write_biom(profile_to_biom(
            data[rank], samples, tree, rankd, named), fp)
    click.echo(' Done.')


def prep_ranks(ranks=None, outmap_dir=None):
    """Prepare directory of output classification maps.

    Parameters
    ----------
    ranks : list, optional
        Target ranks.
    outmap_dir : str, optional
        Path to output file or directory.

    Returns
    -------
    list, dict or None
        Ranks, rank-to-directory map (or None if not necessary).
    """
    # determine ranks
    if ranks is None:
        ranks = ['none']
    else:
        ranks = ranks.split(',')

    # check output directory
    if not outmap_dir:
        return ranks, None
    makedirs(outmap_dir, exist_ok=True)

    # determine output map directory per rank
    if len(ranks) == 1:
        rank2dir = {ranks[0]: outmap_dir}
    else:
        rank2dir = {}
        for rank in ranks:
            dir_ = join(outmap_dir, rank)
            makedirs(dir_, exist_ok=True)
            rank2dir[rank] = dir_
    return ranks, rank2dir
