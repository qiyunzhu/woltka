# Woltka

**Woltka**, or **W**eb **o**f **L**ife **T**ool**k**it **A**pp, is a collection of analytical tools for parsing. It means to interface between shotgun metagenomics (e.g., sequence aligners, taxonomic & functional profilers) and advanced analytical platforms (e.g., QIIME 2).

It currently provides the following functional modules (subcommands):

- **gotu**: gOTU table generation workflow.
- **classify**: Complete classification workflow with all parameters.

Woltka ships with a **QIIME 2 plugin**. [See here for instructions](q2).


## Installation

Requires: Python 3.6 or above.

```bash
pip install git+https://github.com/qiyunzhu/woltka.git@dev
```

After installation, launch the program by executing:

```bash
woltka
```


## Example usage

Woltka ships with small test datasets under this directory:

```
<program_dir>/woltka/tests/data
```

One can execute the following commands to make sure that Woltka functions correctly, and to get an impression of the basic usage of Woltka.

gOTU table generation:

```bash
woltka gotu -i align/bowtie2 -o table.tsv
```

Taxonomic profiling at the ranks of phylum, genus and species as defined in NCBI taxdump:

```bash
woltka classify \
  -i align/bowtie2 \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank phylum,genus,species \
  -o output_dir
```

Functional profiling by UniRef entries then by GO terms (molecular process):

```bash
woltka classify \
  -i align/bowtie2 \
  --coords function/coords.txt.xz \
  --map function/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --rank uniref,process \
  -o output_dir
```


(Contents below are under development.)

## gOTU table generation

## Tree-based classification

Woltka features a highly flexible hierarchical classification system. It is represented by a tree structure instead of a fixed number of levels (e.g., the eight standard taxonomic ranks). In another word, it is **rank-free**.

Wolkta supports various format of classification systems, specifically:

1. `--nodes`: NCBI-style `nodes.dmp` or compatible formats, in which each taxon is pointed to its parent taxon.
2. `--newick`: Newick-format tree, in which labels of tips and internal nodes are considered as taxa.
3. `--ranktb`: Tab-delimited table in which each column represents a rank
4. `--lineage`: Map of taxon to lineage string, in the format of `taxonomic;units;from;high;to;low`. Can be Greengenes-style taxonomy where level codes such as `k__` will be parsed.
  - Compatible with QIIME, SHOGUN, MetaPhlAn2, GTDB, etc.
5. `--map`: One or multiple ordered maps of lower taxa to higher taxa. For example, the 1st file maps genes to UniRef entries, the 2nd maps UniRef entries to GO terms, the 3rd maps GO terms to GO slim terms, so on so forth.
6. If no classification file is provided, Woltka will automatically build a classification system from the mapping files, in which subject identifiers will be parsed as lineage strings.

## Integrating taxonomic & functional analyses

Woltka allows assignment of sequences to both taxonomic units AND functional categories using the same set of alignment files.
