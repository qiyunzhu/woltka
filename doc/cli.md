# Command-line interface

The command-line interface (CLI) of Woltka provides several commands:

- **gotu**: gOTU table generation.
- **classify**: Complete classification workflow with all parameters.

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
`--lines` | Number of lines to read from alignment file per chunk. Default: 1000000.

### Hierarchies

* See [classification system](classify.md) for details.

Option | Description
--- | ---
`--nodes` | Hierarchies defined by NCBI nodes.dmp or compatible formats.
`--newick` | Hierarchies defined by a tree in Newick format.
`--lineage` | Map of lineage strings. Can accept Greengenes-style rank prefix.
`--rank-table` | Table of classification units at each rank (column).
`--map`, `-m` | 'Map(s) of subjects or lower classification units to higher ones. Can accept multiple maps.
`--names` | Names of classification units as defined by NCBI names.dmp or a plain map.

### Assignment

Option | Description
--- | ---
`--rank`, `-r` | Classify sequences at this rank. Ignore or enter "none" to omit classification; enter "free" for free-rank classification. Can specify multiple comma-delimited ranks and one profile will be generated for each rank.
`--above/--no-above` | Allow assigning to a classification unit higher than given rank. Default: auto-decision.
`--major` | Majority-rule assignment percentage threshold. Range: [51, 99].
`--ambig/--no-ambig` | Allow assigning one sequence to multiple classification units. Default: True.
`--subok/--no-subok` | Can report subject IDs in classification result. Default: True.
`--deidx/--no-deidx` | Strip "underscore index" suffixes from subject IDs. Default: False.

### Gene matching

Option | Description
--- | ---
`--coords`, `-c` | Reference gene coordinates on genomes.
`--overlap`, | Read/gene overlapping percentage threshold. Default: 80.

### Output files

Option | Description
--- | ---
`--to-biom/--to-tsv` | Force output feature table format (BIOM or TSV). If omitted, format defaults to BIOM if there are multiple ranks, or based on output filename extension (`.biom` for BIOM, otherwise TSV) if there is one rank.
`--name-as-id` | Replace feature IDs with names.
`--add-rank` | Append feature ranks to table.
`--add-lineage` | Append lineage strings to table.
`--outmap`, `-u` | Write read-to-feature maps to directory.
`--outmap-zip` | Compress read maps using this algorithm. Options: `none`, `gz` (default), `bz2`, `xz`.
