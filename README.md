# Woltka

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Build Status](https://travis-ci.org/qiyunzhu/woltka.svg?branch=master)](https://travis-ci.org/qiyunzhu/woltka)
[![Coverage Status](https://coveralls.io/repos/github/qiyunzhu/woltka/badge.svg?branch=master)](https://coveralls.io/github/qiyunzhu/woltka?branch=master)

**Woltka** (Web of Life Toolkit App), is a bioinformatics package for shotgun metagenome data analysis. It takes full advantage of, and it not limited by, the [WoL](https://biocore.github.io/wol/) reference phylogeny. It bridges first-pass sequence aligners with advanced analytical platforms (such as QIIME 2). Highlights of this program include:

- gOTU: fine-grain community ecology.
- Tree-based, rank-free classification.
- Combined taxonomic & functional analysis.

Woltka ships with a **QIIME 2 plugin**. [See here for instructions](woltka/q2).

## Contents

- [Overview](#overview)
- [Installation](#installation)
- [Example usage](#example-usage)
- Main workflow
  - [Input files](doc/input.md)
  - [Output files](doc/output.md)
  - [Classification systems](doc/hierarchy.md)
  - [Classification methods](doc/classify.md)
  - [Coordinates matching](doc/ordinal.md)
  - [Stratification](doc/stratify.md)
- Profile tools
  - [Collapse](doc/collapse.md), [Coverage](doc/coverage.md), [Filter](doc/filter.md), [Merge](doc/merge.md)
- Tutorials
  - [Working with WoL](doc/wol.md)
  - [gOTU analysis](doc/gotu.md)
- For users of
  - [QIIME 2](woltka/q2), [Qiita](doc/app.md#qiita), [SHOGUN](doc/wol.md#sequence-alignment), [GTDB](doc/gtdb.md), [MetaCyc](doc/metacyc.md)
- References
  - [Command-line interface](doc/cli.md)
  - [Computational efficiency](doc/perform.md)
- [FAQs](#doc/faq.md)
- [Notes](#notes)


## Overview

### Where does Woltka fit in a pipeline

Woltka is a **classifier**. It serves as a middle layer between sequence alignment and community ecology analyses.

### What does Woltka do

Woltka processes **alignments** -- the mappings of query sequences against reference sequences (such as microbial genomes or genes), and infers the best placement of the queries in a hierarchical classification system. One query could have simultaneous matches in multiple references. Woltka finds the most suitable classification unit(s) to describe the query accordingly the criteria specified by the researcher. Woltka generates **profiles** (feature tables) -- the frequencies (counts) of classification units which describe the composition of samples.

### What else does Woltka do

Woltka provides several utilities for handling feature tables, including collapsing a table to higher-level features, calculating feature group coverage, filtering features based on per-sample abundance, and  merging tables.

### What does Woltka not do

Woltka does NOT **align** sequences. You need to align your FastQ (or Fast5, etc.) files against a reference database (we recommend [WoL](https://biocore.github.io/wol/)) use an aligner of your choice (BLAST, Bowtie2, etc.). The resulting alignment files can be fed into Woltka.

Woltka does NOT **analyze** profiles. We recommend using [QIIME 2](https://qiime2.org/) for robust downstream analyses of the profiles to decode the relationships among micobial communities and with their environments.


## Installation

Requirement: [Python](https://www.python.org/) 3.6 or above, with Python package [biom-format](http://biom-format.org/).

```bash
pip install git+https://github.com/qiyunzhu/woltka.git
```

After installation, launch the program by executing:

```bash
woltka
```

[More details about installation](doc/install.md) are provided here.

## Example usage

Woltka provides several small test datasets under [woltka/tests/data](woltka/tests/data). To access them, [download](https://github.com/qiyunzhu/woltka/archive/master.zip) this GitHub repo, unzip, and navigate to this directory.

One can execute the following commands to make sure that Woltka functions correctly, and to get an impression of the basic usage of Woltka.

(Note: a more complete list of commands at provided [here](woltka/tests/data). Alternatively, you can skip this test dataset check out the [instructions](doc/wol.md) for working with WoL.)

1\. gOTU table generation ([details](doc/gotu.md)):

```bash
woltka gotu -i align/bowtie2 -o table.biom
```

The input path, [`align/bowtie2`](woltka/tests/data/align/bowtie2), is a directory containing five Bowtie2 alignment files (`S01.sam.xz`, `S02.sam.xz`,... `S05.sam.xz`) (SAM format, xzipped), each representing the mapping of shotgun metagenomic sequences per sample against a reference genome database.

The output file, `table.biom`, is a feature table in BIOM format, which can then be analyzed using various bioformatics programs such as [QIIME 2](https://qiime2.org/).

2\. Taxonomic profiling at the ranks of phylum, genus and species ([details](doc/hierarchy.md)):

```bash
woltka classify \
  -i align/bowtie2 \
  --map taxonomy/g2tid.txt \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank phylum,genus,species \
  -o output_dir
```

The mapping file (`g2tid.txt`) translates genome IDs to taxonomic IDs, which then allows Woltka to classify query sequences based on the NCBI taxonomy (`nodes.dmp` and `names.dmp`).

The output directory (`output_dir`) will contain three feature tables: `phylum.biom`, `genus.biom` and `species.biom`, each representing a taxonomic profile at one of the three ranks.

3\. Functional profiling by UniRef entries then by GO terms (molecular process):

```bash
woltka classify \
  -i align/bowtie2 \
  --coords function/coords.txt.xz \
  --map function/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --map-as-rank \
  --rank uniref,process \
  -o output_dir
```

Here, the input files are still read-to-genome alignments, instead of read-to-gene ones, but Woltka matches reads to genes based on their coordinates on the genomes (as indicated by the file `coords.txt`). This ensures consistency between taxonomic and functional classifications.

Subsequently, Woltka is able to assign query sequences to functional units, as defined in mapping files (`uniref.map` and `go/process.tsv`). As you can see, compressed files are supported and auto-detected.

Similarly, the output files are two functional profiles: `uniref.biom` and `process.biom`.

One can also combine taxonomic and functional profilings in a **stratification** analysis. See [details](doc/stratify.md).

## Notes

### Citation

Woltka is currently under development. Please directly cite this GitHub repository:

- https://github.com/qiyunzhu/woltka

### Grants

The development of Woltka is supported by: (to be added).

### Contact

Please forward any questions to the project leader: Dr. Qiyun Zhu (qiz173@ucsd.edu) or the senior PI: Dr. Rob Knight (robknight@ucsd.edu).
