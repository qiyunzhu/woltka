# Tree-based classification

Woltka features a highly flexible hierarchical classification system. It is represented by a **tree** structure instead of a fixed number of levels (e.g., the eight standard taxonomic ranks). In another word, it is **rank-free**.

Wolkta supports various formats of classification systems, specifically:

1. `--nodes`: NCBI-style `nodes.dmp` (columns are delimited by "\<tab\>|\<tab\>") or a tab-delimited file, in which each taxon (1st column) points to its parent taxon (2nd column). Rank (3rd column) is optional.

2. `--newick`: Newick-format tree, in which labels of nodes (tips, internal nodes and root) are considered as taxa. All nodes must have labels and all labels must be unique.

3. `--ranktb`: Table of per-taxon per-rank assignments. Each column represents a rank. Column header will be treated as rank name.

4. `--lineage`: Map of taxon to lineage string (`;`-delimited taxa from high to low)

   It can be **Greengenes**-style taxonomy where rank codes such as `k__` will be parsed. But the rank code is not mandatory. Unassigned taxon (e.g., `s__`) and non-unique taxon are acceptable (e.g., `p__Actinobacteria` and `c__Actinobacteria`).

   Compatible with widely-used taxonomy systems in e.g., QIIME, SHOGUN, MetaPhlAn2, GTDB, etc.

5. `--map` or `-m`: Simple map of lower taxon \<tab\> higher taxon.

   One can supply multiple maps (by entering multiple `--map` parameters) to constitute several hierarchies. For example, the 1st file maps genes to UniRef entries, the 2nd maps UniRef entries to GO terms, the 3rd maps GO terms to GO slim terms, so on so forth.

   Flag `--map-is-rank` is to instruct the program to treat the map filename as rank. For example, with this flag, taxa in the 2nd column of `uniref.map.gz` will be given the rank "uniref".

If no classification file is provided, Woltka will automatically build a classification system from the alignment files, in which subject identifiers will be parsed as lineage strings.

Subjects themselves are part of the classification system. A map of subject to taxon (e.g., a genome ID to NCBI TaxID map) can be supplied with the `--map` parameter if necessary.

Classification files are **additive**, i.e., if multiple files of the same or different formats are provided, all of them will be parsed and added to the classification hierarchies -- unless they conflict -- which will be noted by Woltka. Example command (example files provided under [`taxonomy`](woltka/tests/data/taxonomy), same below):

```bash
woltka classify \
  --map nucl2g.txt \
  --map g2taxid.txt \
  --nodes taxdump/nodes.dmp \
  --names taxdump/names.dmp \
  ...
```

In this command, three layers of hierarchies are provided: 1) nucleotide ID to genome ID (`nucl2g.txt`), 2) genome ID to taxonomy ID (`g2taxid.txt`), 3) NCBI taxonomy tree (`nodes.dmp`).

Again, compressed files are supported and automatically recognized. The following command works:

```bash
woltka classify \
  --lineage gg_13_5_taxonomy.txt.gz
  ...
```

Furthermore, one can supply Woltka with a taxon name dictionary, and the output profile will show taxon names instead of taxon IDs:

* `--names`: NCBI-style `names.dmp` or a simple map of taxon \<tab\> name. Example:

```bash
woltka classify \
  --nodes taxdump/nodes.dmp \
  --names taxdump/names.dmp \
  ...
```

## Combined taxonomic & functional analyses

Woltka combines the two fundamental analyses in metagenomics: taxonomic profiling (mapping reads to genomes) and functional profiling (mapping reads to functional genes) into one run. This saves compute, ensures consistency, and allows for stratification which is essential for understanding functional diversity across the microbial community.

This is achieved by an efficient algorithm implemented in Woltka, which matches read alignments and annotated genes based on their coordinates on the host genome.

The coordinates of read-to-genome alignments are provided in the alignment files. One needs to provide Woltka with a table of gene coordinates. The format is like:

WoL reannotations (available for [download](https://biocore.github.io/wol/)):

```
>G000006745
1       806     372
2       2177    816
3       3896    2271
4       4446    4123
5       4629    4492
```

Native NCBI annotations and accessions:

```
## GCF_000005825.2
# NC_013791.2
WP_012957018.1  816 2168
WP_012957019.1  2348    3490
WP_012957020.1  3744    3959
WP_012957021.1  3971    5086
...
```

Again, compressed files are supported.

With the coordinates file, one can streamline the read-to-gene matching step into a Woltka protocol. Here is an example for functional profiling:

```bash
woltka classify \
  --coords coords.txt \
  --map gene2function.txt \
  --map function2pathway.txt \
  ...
```
