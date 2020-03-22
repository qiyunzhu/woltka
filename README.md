# Woltka

[![Build Status](https://travis-ci.org/qiyunzhu/woltka.svg?branch=master)](https://travis-ci.org/qiyunzhu/woltka)
[![Coverage Status](https://coveralls.io/repos/github/qiyunzhu/woltka/badge.svg?branch=master)](https://coveralls.io/github/qiyunzhu/woltka?branch=master)

**Woltka** (Web of Life Toolkit App), is a bioinformatics package for shotgun metagenome data analysis. It takes full advantage of, and it not limited by, the [WoL](https://biocore.github.io/wol/) reference phylogeny. It bridges first-pass sequence aligners with advanced analytical platforms (such as QIIME 2). Highlights of this program include:

- gOTU: fine-grain community ecology.
- Tree-based, rank-free classification.
- Combined taxonomic & functional analysis.

Woltka ships with a **QIIME 2 plugin**. [See here for instructions](woltka/q2).

## Contents

- [Installation](#installation)
- [Example usage](#example-usage)
- Details
  - [Input files](doc/input.md)
  - [Output files](doc/output.md)
  - [Classification system](doc/classify.md)
  - [Ordinal matching](doc/ordinal.md)
  - [Stratification](doc/stratify.md)
- Tutorials
  - [gOTU analysis](doc/gotu.md)
  - Tree-based classification
  - Combined taxonomic & functional analyses
- For users of
  - [QIIME 2](woltka/q2), [Qiita](doc/app.md#qiita), SHOGUN
  - Bowtie2, BWA, Minimap2
  - BLAST, DIAMOND, VSEARCH
- References
  - [Command-line interface](doc/cli.md)
  - [Parameter auto-decision](doc/auto.md)
  - [Developer's guidlines](doc/dev.md)
- [Notes](#notes)

## Installation

Requirement: Python 3.6 or above.

```bash
pip install git+https://github.com/qiyunzhu/woltka.git
```

After installation, launch the program by executing:

```bash
woltka
```

## Example usage

Woltka provides several small test datasets under [woltka/tests/data](woltka/tests/data). To access them, [download](https://github.com/qiyunzhu/woltka/archive/master.zip) this GitHub repo, unzip, and navigate to this directory.

One can execute the following commands to make sure that Woltka functions correctly, and to get an impression of the basic usage of Woltka.

1\. gOTU table generation ([details](doc/gotu)):

```bash
woltka gotu -i align/bowtie2 -o table.biom
```

The input path, [`align/bowtie2`](woltka/tests/data/align/bowtie2), is a directory containing five Bowtie2 alignment files (`S01.sam.xz`, `S02.sam.xz`,... `S05.sam.xz`) (SAM format, xzipped), each representing the mapping of shotgun metagenomic sequences per sample against a reference genome database.

The output file, `table.biom`, is a feature table in BIOM format, which can then be analyzed using various bioformatics programs such as [QIIME 2](https://qiime2.org/).

2\. Taxonomic profiling at the ranks of phylum, genus and species ([details](doc/classify)):

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

## Notes

### Citation

Woltka is currently under development. Please directly cite this GitHub repository:

- https://github.com/qiyunzhu/woltka

### Grants

The development of Woltka is supported by: (to be added).

### Contact

Please forward any questions to the project leader: Dr. Qiyun Zhu (qiz173@ucsd.edu) or the senior PI: Dr. Rob Knight (robknight@ucsd.edu).
