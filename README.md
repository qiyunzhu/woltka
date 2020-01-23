# Woltka

**Woltka**, or **W**eb **o**f **L**ife **T**ool**k**it **A**pp, is a collection of analytical tools for parsing. It means to interface between shotgun metagenomics (e.g., sequence aligners, taxonomic & functional profilers) and advanced analytical platforms (e.g., QIIME 2).


## gOTU table generation


## Tree-based classification

Woltka features a highly flexible hierarchical classification system. It is represented by a tree structure instead of a fixed number of levels (e.g., the eight standard taxonomic ranks). In another word, it is **rank-free**.

Wolkta supports various format of classification systems, specifically:

1. `--nodes`: NCBI-style `nodes.dmp` or compatible formats, in which each taxon is pointed to its parent taxon.
2. `--newick`: Newick-format tree, in which labels of tips and internal nodes are considered as taxa.
3. `--levels`: Tab-delimited table in which each column represents a level
4. `--lineage`: Map of taxon to lineage string, in the format of `taxonomic;units;from;high;to;low`. Can be Greengenes-style taxonomy where level codes such as `k__` will be parsed.
5. `--groups`: One or multiple ordered maps of lower taxa to higher taxa. For example, the 1st file maps genes to UniRef entries, the 2nd maps UniRef entries to GO terms, the 3rd maps GO terms to GO slim terms, so on so forth.
6. If no classification file is provided, Woltka will automatically build a classification system from the mapping files, in which subject identifiers will be parsed as lineage strings.


## Integrating taxonomic & functional analyses

Woltka allows assignment of sequences to both taxonomic units AND functional categories using the same set of alignment files.