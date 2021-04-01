# Working with "Web of Life" (WoL)

The "Web of Life" (WoL) project is a series efforts to reconstruct an accurate reference phylogeny for microbial genomes, and to build resources that can (and are already) benefiting microbiome researchers.

Phase I of the project was already completed ([Zhu et al., 2019](https://www.nature.com/articles/s41467-019-13443-4)). We have released a reference [tree](https://biocore.github.io/wol/data/trees/astral/branch_length/cons/collapsed/astral.cons.nid.e5p50.nwk), built on 10,575 bacterial and archaeal genomes, based on 381 marker genes.

The project is detailed at our website: https://biocore.github.io/wol/, including data and [metadata](https://biocore.github.io/wol/data/genomes/metadata.tsv.bz2), code, protocols, a gallery and a visualizer. Large data files are hosted at our Globus endpoint: [WebOfLife](https://app.globus.org/file-manager/collections/31acbeb8-c62f-11ea-bef9-0e716405a293) (see [instruction](https://biocore.github.io/wol/download#download-via-globus)).

This public resource provides everything one needs to start microbiome data analysis using WoL, including raw sequence data, metadata, tree and taxonomy, and pre-built databases that are ready to be plugged into your bioinformatics protocols. Currently, we provide databases for QIIME 2, SHOGUN, Bowtie2, Centrifuge, Kraken2 / Bracken, BLASTn and BLASTp, Minimap2, and DIAMOND. Even if your favorate tool is not on this list, we provide detailed tutorials on how to [build your own database](https://biocore.github.io/wol/protocols/genome_database) and many other related [protocols](https://biocore.github.io/wol/protocols/). Meanwhile, WoL is also hosted at our web-based microbiome study platform: **Qiita** (https://qiita.ucsd.edu/) (see [details](qiita.md)).

The following tutorial assume that you have downloaded the entire WoL directory from our Globus server. The paths mentioned below are relative to this directory.

## Sequence alignment

First, you need to align your sequences (namely your FastQ / Fast5 / BAM files) against the WoL database using an aligner of your choice. Let's take [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for example. Our bioinformatics tool, [SHOGUN](https://github.com/knights-lab/SHOGUN), provides a Bowtie2 wrapper optimized for shotgun metagenomic datasets:

```bash
shogun align -d databases/shogun -a bowtie2 -t 16 -p 0.95 -i input.fa -o .
```

This will generate a [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) format alignment file.

The alignment step has been automated in [Qiita](https://qiita.ucsd.edu/). If you use Qiita, the SAM file is ready for download.

[**Note**] You can also run Bowtie2 manually using your choice of parameters, or using other aligners and other databases. Woltka is designed for flexibility.

## OGU analysis

```bash
woltka classify -i input.sam -o output.biom
```

Note that you can [compress](perform.md#compress-alignment-files) the SAM file to save disk space, and Woltka can parse compressed files.


## Free-rank classification

Use the original NCBI taxonomy:

```bash
woltka classify \
  --input input.sam \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --output output.biom
```

Use lineage strings extracted from NCBI (will lose some resolution, but results are more structured, especially for users familiar with QIIME 2):

```bash
woltka classify \
  --input input.sam \
  --lineage taxonomy/lineage.txt \
  --output output.biom
```

We also provide original and curated [NCBI](https://biocore.github.io/wol/data/taxonomy/ncbi/) and [GTDB](https://biocore.github.io/wol/data/taxonomy/gtdb/) taxonomy for choice.


## Classification at specific ranks

Slightly modify the command, adding desired ranks:

```bash
woltka classify \
  --input input.sam \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank phylum,genus,species \
  --output output.biom
```

## Coordinates-based functional classification

Associate read alignments with ORFs, and move up multiple functional levels (protein, reaction, pathway...) using a cascade of mapping files. The example below uses the MetaCyc system ([see details](metacyc.md)):

```bash
mcdir=annotation/metacyc
woltka classify \
  --input input.sam \
  --coords annotation/coords.txt.xz \
  --map $mcdir/wol-to-protein.txt.xz --names $mcdir/protein_name.txt \
  --map $mcdir/protein-to-enzrxn.txt --names $mcdir/enzrxn_name.txt \
  --map $mcdir/enzrxn-to-reaction.txt --names $mcdir/reaction_name.txt \
  --map $mcdir/reaction-to-pathway.txt --names $mcdir/pathway_name.txt \
  --map-as-rank \
  --rank protein,enzrxn,reaction,pathway \
  --output output_dir
```

Note: This command won't handle multiple mapping (e.g., one protein involved in three pathways). A capable solution is provided [here](metacyc.md).

## Stratified taxonomic / functional classification

Say, you want to stratify functional annotations by genus (taxonomy). First, run taxonomic classification at the genus level, and export read-to-genus maps:

```bash
woltka classify \
  --input input.sam \
  ...
  --rank genus \
  --name-as-id \
  --output genus.biom
  --outmap map_dir
```

Second, run functional annotation, adding the read-to-genus maps for stratification:

```bash
woltka classify \
  --input input.sam \
  --coords annotation/coords.txt.xz \
  ...
  --stratify map_dir
  --output output_dir
```
