# Output files

The main output file(s) of Woltka are **feature tables**, in which column headers are sample IDs and row headers are feature IDs, and each cell value represents the _count_ of a particular feature in a particular sample.

The term "**feature**" can mean various things in different applications. In gOTU analysis, it is a gOTU. In a taxonomic classification analysis, they are taxonomic units at designated ranks, or in rank-free style. In a functional classification analysis, they are


sequence **alignment** files. The term "alignment" here describes the operation of aligning short sequencing reads against long reference sequences. The basic information an alignment file provides is the **mapping** between queries and subjects. In addition, the position and quality of alignments are available in some formats, which are useful in some applications.

## Contents

- [BIOM format](#biom-format)
- [TSV format](#tsv-format)
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
