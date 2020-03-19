# Input files

The input files for Woltka are sequence **alignment** files. The term "alignment" here describes the operation of aligning short sequencing reads against long reference sequences. The basic information an alignment file provides is the **mapping** between queries and subjects. In addition, the position and quality of alignments are available in some formats, which are useful in some applications.

## Contents

- [File formats](#file-formats)
- [Filename patterns](#filename-patterns)
- [Demultiplexing](#Demultiplexing)

## File formats

Woltka supports the following alignment formats (specified by parameter `--format` or `-f`):

- `map`: Simple mapping of query \<tab\> subject.
- `sam`: SAM format. Supported by multiple tools such as Bowtie2 and BWA.
- `b6o`: BLAST tabular format (i.e., BLAST parameter `-outfmt 6`). Supported by multiple tools such as BLAST, DIAMOND, VSEARCH, BURST, etc.

If not specified, Woltka will automatically infer the format of input alignment files.

Woltka supports and automatically detects common file compression formats including `gzip`, `bzip2` and `xz`.

## Filename patterns

In the default mode, Woltka treats every file under the directory specified by the `--input` or `-i` parameter as an alignment file for one sample. The sample ID is automatically extracted from the filename, with filename extension stripped. Compression file extensions are automatically recognized. For example, the sample IDs for `S01.sam` and `S02.m8.gz` are `S01` and `S02`, respectively.

One may restrict this behavior to avoid confusions (e.g., there are unwanted files) by the following two parameters:

- `--filext` or `-e`: A suffix to be stripped from filenames, and the remaining part is considered a sample ID.

  For example, if valid alignment filenames have the pattern of `ID_L001.aln.bz2`, one may specify `-e _L001.aln.bz2`, so that only `ID` is retained, and filenames which do not have this suffix (e.g., `ID.bt2.log` or`readme.txt`) are ignored.

- `--sample-ids` or `-s`: A file containing valid sample IDs (one ID per line) to include in the current analysis. This list also specifies the order of sample IDs in the output profile.

  Only the first column before \<tab\> is considered. Lines starting with `#` are omitted. Therefore, a metadata table may also serve as a valid sample ID list.

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
  --sample-ids ids.txt \
  --no-demux \
  --output profile.biom \
  ...
```
