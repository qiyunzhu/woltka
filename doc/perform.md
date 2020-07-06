# Computational efficiency

Woltka is for handling very large alignment files and complex classification systems. Therefore computational expense is an important consideration. This page explains how much time does Woltka spend on typical datasets and how much memory does it need, as well as multiple tips and tricks for improving performance.


## Contents

- [Estimate runtime and memory](#estimate-runtime-and-memory)
- [Benchmarks on a typical dataset](#benchmarks-on-a-typical-dataset)
- [Compress alignment files](#compress-alignment-files)
- [Convert alignments to simple maps](#convert-alignments-to-simple-maps)


## Estimate runtime and memory

Rules of thumb:

- Woltka's runtime is **linear** (_O_(_n_)) vs. the size of the input alignments, plus a small, constant overhead of reading the database.
- Woltka's memory consumption is **constant**. It is only determined by the size of the database, but is independent from the input alignments.
- Woltka analyses are I/O-intensive. Performance of disk and memory matter.
- Coordinates-based classification takes several folds more time than simple classification.
- Woltka is a **single**-threaded program at the moment.

For the Web of Life ([WoL](https://biocore.github.io/wol/)) database, which includes 10,575 genomes, the memory consumption is **~120 MB** for taxonomic classification, and **~22 GB** for gene coordinates-based functional classification.

## Benchmarks on a typical dataset

In this example, we started with the [CAMI](https://data.cami-challenge.org/) high complexity toy dataset, which contains 5 samples with 15 Gbp HiSeq sequencing data each. We aligned them using the [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) aligner through [SHOGUN](https://github.com/knights-lab/SHOGUN), against the Web of Life ([WoL](https://biocore.github.io/wol/)) database. This step generates up to 16 high-score alignments (matches) for each query sequence.

The five resulting alignment files (SAM format, gzipped) total **38.6 GB** in size. They include **1.05 billion** alignments (lines), with 501 million unique query sequences. They were stored in a hard disk drive (59,000 rpm, SATA 3). The test was performed on with a Xeon E3-1230 v3 CPU (Haswell microarchitecture) and 32 GB DDR3 memory. The software environment was Ubuntu 18.04 plus Python 3.7.

Here are the benchmarks of multiple typical Woltka analyses:

Analysis | Runtime (h:m:s) | Memory
--- | --- | ---
gOTU (no classification) | 1:17:43 | 102.7 MB
Free-rank taxonomic classification | 1:27:43 | 108.0 MB
Taxonomic classification at genus | 1:20:32 | 107.5 MB
Taxonomic classification at genus, while writing gzipped read maps | 1:33:21 | 107.8 MB
Taxonomic classification at three ranks: phylum, genus and species | 1:46:25 | 118.0 MB
Coordinates-based functional classification | 5:36:49 | 20.40 GB


## Compress alignment files

The input files for Woltka -- sequence alignment files cost disk space, and compression can help. An advantage of Woltka is that it automatically detects and parses compressed files, including the three most common formats: **gzip**, **bzip2** and **xz** (i.e., LZMA) (see [details](input.md#file-formats)).

Below is a set of simple benchmarks performed on a typical shotgun metagenomics dataset, including the original FastQ sequence files and two alignment file formats.

![Compression method benchmarks](img/zip.png)

It is evident that most rigorous method xz provides highest compression ratio, though with much longer runtime when compressing and _modestly_ longer runtime when decompressing, compared to the fastest method gzip. If disk space is your main concern and CPU hour is second (remember: compression is a one-off operation), we would recommend adopting **xz** as a rule of thumb.


## Convert alignments to simple maps

Woltka allows you to generate maps of query sequences to classification units via the `--outmap` or `-u` parameter (see [details](output.md#output-read-maps)). These maps themselves can serve as input files for Woltka, enabling efficient reuse AND the powerful [stratification](stratify.md) function.

It is okay to convert alignment files into maps, and ONLY keep the maps while deleting the original alignments from your valuable storage space (but see warning below). Here is how:

```
woltka classify -i input_dir -o whatever.biom -u map_dir
```

The outcome will be per-sample maps of query sequences to original subject sequences (i.e., no classification is involved). In the future, you can repeatedly do this, with choices of parameter sets:

```
woltka classify -i map_dir ...
```

[**Temporary notice**] At this moment Woltka cannot take non-unique mappings (one query to multiple subjects) as input. Therefore, if your alignment files may contain non-unique matches, you need to add `--uniq` to the first command.

Woltka supports outputing compressed maps. The default compression method is **gzip**, and it can be customized via the `--zipmap` paramter (see [details](output.md#output-read-maps)).

The following benchmarks are based on five already-**xz**'ed Bowtie2 SAM files, totalling 327.224 MB:

Method | Runtime (min) | Map size (MB)
--- | --- | ---
Don't generate maps | 1:15.94 | -
Don't compress maps | 1:22.45 | 159.832
Using **gzip** | 1:38.42 | 17.960
Using **bzip2** | 1:39.99 | 14.496
Using **xz** | 2:33.65 | 10.664

It is your call how to balance runtime and map file size.

Warning however, converting alignments to maps will **lose alignment coordinates**, and **blocks coordinates-based analyses** (such as functional classification) (see details). If these analyses are in your plan, you should keep the alignment files.
