# Tips for high performance computing

## Contents

- [Compress alignment files](#compress-alignment-files)
- [Convert alignments to simple maps](#convert-alignments-to-simple-maps)


## Compress alignment files

The iput files for Woltka -- sequence alignment files cost disk space, and compression can help. An advantage of Woltka is that it automatically detects and parses compressed files, including the three most common formats: **gzip**, **bzip2** and **xz** (i.e., LZMA) (see [details](input.md#file-formats)).

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

[**Temporary notice**] At this moment Woltka cannot take non-unique mappings (one query to multiple subjects) as input. Therefore, if your alignment files may contain non-unique matches, you need to add `--no-ambig` to the first command.

Woltka supports outputing compressed maps. The default compression method is **gzip**, and it can be customized via the `--outmap-zip` paramter (see [details](output.md#output-read-maps)).

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
