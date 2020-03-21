# Input files

The input files for Woltka are sequence **alignment** files. The term "alignment" here describes the operation of aligning short sequencing reads against long reference sequences. The basic information an alignment file provides is the **mapping** between queries and subjects. In addition, the position and quality of alignments are available in some formats, which are useful in some applications.

## Contents

- [Input filepath](#input-filepath)
- [File formats](#file-formats)
- [Filename pattern](#filename-pattern)
- [Sample list](#sample-list)
- [Demultiplexing](#demultiplexing)

## Input filepath

Parameter `--input` or `-i` is to let Woltka know where to find the input alignment file(s). It can be any of the following three scenarios:

1\. A **directory**. Woltka will search this directory for alignment files, and treat each of them as _one sample_. One may control. See below for [details](#filename-patterns).

2\. A single **alignment file**. Woltka considers it as a _multiplexed_ alignment file. Sample IDs will be extracted from sequence identifiers (see [demultiplexing](#demultiplexing) below).

3\. A **mapping file** of sample ID \<tab\> alignment file path. The paths must point to existing files. They can either be full paths, or simple filenames under the same directory as the mapping file. For example, one can place a `map.txt` of the following content to where alignment files are located.

```
S01 S01.sam
S02 S02.sam
S03 S03.sam
...
```

## File formats

Woltka supports the following alignment formats (specified by parameter `--format` or `-f`):

- `map`: A **simple map** in the format of query \<tab\> subject.
- `sam`: **SAM** format. Supported by multiple tools such as Bowtie2 and BWA.
- `b6o`: **BLAST** tabular format (i.e., BLAST parameter `-outfmt 6`). Supported by multiple tools such as BLAST, DIAMOND, VSEARCH, BURST, etc.

  If not specified, Woltka will automatically infer the format of input alignment files.

Woltka supports and automatically detects common file compression formats including `gzip`, `bzip2` and `xz`.

## Filename pattern

When Woltka searches a directory specified by `--input`, by default, it considers every file as an alignment file for one sample. The sample ID is automatically extracted from the filename, with filename extension stripped. Compression file extensions are automatically recognized. For example, the sample IDs for `S01.sam` and `S02.m8.gz` are `S01` and `S02`, respectively.

One may restrict this behavior to avoid confusions. One method is to use parameter `--filext` or `-e`. It is a suffix to be stripped from filenames, and the remaining part is considered a sample ID.

  For example, if valid alignment filenames have the pattern of `ID_L001.aln.bz2`, one may specify `-e _L001.aln.bz2`, so that only `ID` is retained, and filenames which do not have this suffix (e.g., `ID.bt2.log` or`readme.txt`) are ignored.

A second method is detailed below.

## Sample list

Parameter `--samples` or `-s` is to instruct Woltka which samples are to be included in the analysis, and in which order they should appear in the output feature table. If not specified, all samples will appear in alphabetical order. This feature applies to all three types of `--input` parameters (see [above](#input-filepath)).

  It can be a list of sample IDs separated by comma (e.g., `S01,S02,S03...`), if there aren't many IDs to type. Or,

  It can point to a file containing sample IDs, one ID per line. Only the first column before \<tab\> is considered. Lines starting with `#` are omitted. Therefore, a metadata table may also serve as a valid sample ID list. For example:

```
#SampleID
S01
S02
S03
...
```

## Demultiplexing

Woltka supports convenient demultiplexing. If the input path (specified by `--input` or `-i`) points to a file instead of a directory, Woltka will treat it as a multiplexed alignment file.

Specifically, the program divides sample ID and read ID by the first underscore in each query identifier in the alignment file (i.e., a pattern of `sampleID_readID`).

One may manually switch on or off the demultiplexing function by adding flag `--demux` or `--no-demux` to the command. This is useful when there are several multiplexed alignment files (e.g., each from one sequencing lane) in one **directory**, or the only input alignment **file** provided is not multiplexed but just for a single sample.

Example of a complete command:

```bash
woltka classify \
  --input blast_output/ \
  --format b6o \
  --filext .blast6out.gz \
  --samples ids.txt \
  --no-demux \
  --output profile.biom \
  ...
```
