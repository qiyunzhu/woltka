# Command-line interface

The command-line interface (CLI) of Woltka provides several commands:

- [**gotu**](#gotu): gOTU table generation.
- [**classify**](#classify): Complete classification workflow with all parameters.
- [**tools**](#tools): Utilities for working with alignments, maps and profiles.

## Classify

### Basic

Option | Description
--- | ---
`--input`, `-i` (required) | Path to input alignment file or directory of alignment files.
`--output`, `-o` (required) | Path to output profile file or directory of profile files.

### Input files

* See [input files](input.md) for details.

Option | Description
--- | ---
`--format`, `-f` | Format of read alignments. Options: <ul><li>`b6o`: [BLAST tabular format](https://www.ncbi.nlm.nih.gov/books/NBK279684/).</li><li>`sam`: [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).</li><li>`map`: A simple map of query \<tab\> subject</li></ul>If not specified, program will automatically infer from file content.
`--filext`, `-e` | Input filename extension following sample ID.
`--samples`, `-s` | List of sample IDs to be included.
`--demux/--no-demux` | Demultiplex alignment by first underscore in query identifier.
`--trim-sub` | Trim subject IDs at the last given delimiter.

### Hierarchies

* See [classification system](hierarchy.md) for details.

Option | Description
--- | ---
`--nodes` | Hierarchies defined by NCBI nodes.dmp or compatible formats.
`--newick` | Hierarchies defined by a tree in Newick format.
`--lineage` | Lineage strings. Can accept Greengenes-style rank prefix.
`--columns` | Table of classification units per rank (column).
`--map`, `-m` | Mapping of lower classification units to higher ones.
`--map-as-rank` | Extract rank name from mapping filename.
`--names`, `-n` | Names of classification units as defined by NCBI names.dmp or a plain map.

### Assignment

Option | Description
--- | ---
`--rank`, `-r` | Classify sequences at this rank. Enter "none" to directly report subjects; enter "free" for free-rank classification.; enter "free" for free-rank classification. Can specify multiple comma-delimited ranks and one profile will be generated for each rank. If omitted, the program will do "free" if a classification system is provided or "none" if not.
`--uniq` | One sequence can only be assigned to one classification unit, or remain unassigned if there is ambiguity. Otherwise, all candidate units are reported and their counts are normalized.
`--major` | In given-rank classification, use majority rule at this percentage threshold to determine assignment when there are multiple candidates. Range: [51, 99]. Overrides "above" and "uniq".
`--above` | In given-rank classification, allow assigning a sequence to a higher rank if it cannot be assigned to the current rank. Overrides "uniq".
`--subok` | In free-rank classification, allow assigning a sequence to its direct subject, if applicable, before going up in hierarchy.

### Gene matching

Option | Description
--- | ---
`--coords`, `-c` | Reference gene coordinates on genomes.
`--overlap`, | Read/gene overlapping percentage threshold. Default: 80.

### Stratification

Option | Description
--- | ---
`--stratify`, `-t` | Directory of read-to-feature maps for stratification. One file per sample.

### Output files

Option | Description
--- | ---
`--to-biom/--to-tsv` | Force output profile format (BIOM or TSV). If omitted, format defaults to BIOM if there are multiple ranks, or based on output filename extension (`.biom` for BIOM, otherwise TSV) if there is only one rank.
`--unassigned` | Report unassigned sequences (will be marked as "Unassigned").
`--name-as-id` | Replace feature IDs with names. Otherwise append names to table as a metadata column.
`--add-rank` | Append feature ranks to table as a metadata column.
`--add-lineage` | Append lineage strings to table as a metadata column.
`--outmap`, `-u` | Write read-to-feature maps to this directory.
`--zipmap` | Compress read-to-feature maps using this algorithm. Options: `none`, `gz` (default), `bz2`, `xz`.

### Performance

Option | Description
--- | ---
`--chunk` | Number of alignment lines to read and parse in each chunk. Default: 1,000 for plain mapping, or 1,000,000 for ordinal mapping.
`--cache` | Number of recent classification results to cache for faster subsequent classifications. Default: 1024.
`--no-exe` | Disable calling external programs (`gzip`, `bzip2` and `xz`) for decompression. Otherwise, Woltka will use them if available for faster processing, or switch back to Python if not.


## Tools

### Filter

Filter a profile by **per-sample** abundance.

Option | Description
--- | ---
`--input`, `-i` (required) | Path to input alignment file or directory of alignment files.
`--output`, `-o` (required) | Path to output profile file or directory of profile files.
`--min-count`, `-c` | Per-sample minimum count threshold (>=1).
`--min-percent`, `-p` | Per-sample minimum percentage threshold (<100).

### Merge

Merge multiple profiles into one profile.

Option | Description
--- | ---
`--input`, `-i` (required) | Path to input profiles or directories containing profiles. Can accept multiple paths.
`--output`, `-o` (required) | Path to output profile.

### Collapse

Collapse a profile based on feature mapping (supports **many-to-many** mapping).

Option | Description
--- | ---
`--input`, `-i` (required) | Path to input profile.
`--map`, `-m` (required) | Path to mapping of source features to target features.
`--output`, `-o` (required) | Path to output profile.
`--normalize`, `-z` | Count each target feature as 1 / _k_ (_k_ is the number of targets mapped to a source). Otherwise, count as one.
`--names`, `-n` | Path to mapping of target features to names. The names will be appended to the collapsed profile as a metadata column.
