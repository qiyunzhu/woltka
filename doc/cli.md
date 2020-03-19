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
`--outmap` | Write per-sample per-read classification maps to this directory.

### Input files

* See [input files](input) for details.

Option | Description
--- | ---
`--format`, `-f` | Format of read alignments. Options: <ul><li>`b6o`: BLAST tabular format.</li><li>`sam`: SAM format.</li><li>`map`: Simple map of query \<tab\> subject</li></ul>If not specified, program will automatically infer from file content.
`--filext`, `-e` | Input filename extension following sample ID.
`--sample-ids`, `-s` | List of sample IDs to be included.
`--demux/--no-demux` | Demultiplex alignment by first underscore in query identifier.

### classification

* See [classification system](classify) for details.

Option | Description
--- | ---
`--rank`, `-r` | Classify sequences at this rank. Ignore or enter "none" to omit classification; enter "free" for free-rank classification. Can specify multiple comma-delimited ranks and one profile will be generated for each rank.
`--above/--no-above` | Allow assigning to a classification unit higher than given rank. Default: auto-decision.
`--major` | Majority-rule assignment percentage threshold. Range: [51, 99].
`--ambig/--no-ambig` | Allow assigning one sequence to multiple classification units. Default: True.
`--subok/--no-subok` | Can report subject IDs in classification result. Default: True.
`--deidx/--no-deidx` | Strip "underscore index" suffixes from subject IDs. Default: False.

### gene information

Option | Description
--- | ---
`--coords`, `-c` | Table of gene coordinates of  on reference genomes.
`--overlap`, | Read/gene overlapping percentage threshold. Default: 80.

### tree information

Option | Description
--- | ---
`--names` | Names of classification units as defined by NCBI names.dmp or a plain map.
`--nodes` | Hierarchies defined by NCBI nodes.dmp or compatible formats.
`--newick` | Hierarchies defined by a tree in Newick format.
`--lineage` | Map of lineage strings. Can accept Greengenes-style rank prefix.
`--rank-table` | Table of classification units at each rank (column).
`--map`, `-m` | 'Map(s) of subjects or lower classification units to higher ones. Can accept multiple maps.
`--map-is-rank` | Map filename stem is rank name.

### performance
`--lines` | Number of lines to read from alignment file per chunk. Default: 1000000.
