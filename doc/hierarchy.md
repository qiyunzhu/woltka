# Tree-structured classification system

Woltka features a highly flexible hierarchical classification system. It is represented by a **tree** structure, instead of a fixed number of levels (e.g., the eight standard [taxonomic ranks](https://en.wikipedia.org/wiki/Taxonomic_rank)). It can be the traditional taxonomy, or an actual phylogenetic tree, or the hierarchies of gene-module-pathway, or any hierarchical representations of feature relationships.

The term "**rank**" (or "level") is still relevant, but it is merely a property of a feature, and does not bear information of hierarchy. For example, above a _genus_-level unit it does not have to be a _family_-level one, but could directly go to _order_, or have a _tribe_ which isn't common for the rest of the tree, or one or more nodes which do NOT have the rank assignment.

In another word, Woltka classification is **rank-free**. This design enables finer-grain resolution of feature relationships, in addition to flexibility. It is therefore suitable for complex systems, such as phylogenetic trees.

That being said, Woltka still supports ranked hierarchies and one can instruct the program to target one or more specific ranks.


## Contents

- [Supported hierarchy files](#supported-hierarchy-files)
- [How Woltka handles hierarchies](#how-woltka-handles-hierarchies)
- [Feature name dictionary](#feature-name-dictionary)


## Supported hierarchy files

Wolkta supports various types and formats of classification systems, specifically:

1. `--nodes`: [NCBI-style](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) `nodes.dmp` (columns are delimited by "\<tab\>|\<tab\>") or a tab-delimited file, in which each taxon (1st column) points to its parent taxon (2nd column). Rank (3rd column) is optional.

2. `--newick`: [Newick-format](https://en.wikipedia.org/wiki/Newick_format) tree, in which labels of nodes (tips, internal nodes and root) are considered as taxa. All nodes must have labels and all labels must be unique.

3. `--columns`: Table of per-taxon per-rank assignments. Each column represents a rank. Column header will be treated as rank name.

4. `--lineage`: Map of taxon to lineage string (`;`-delimited taxa from high to low)

   It can be [Greengenes-style](http://greengenes.secondgenome.com/) taxonomy where rank codes such as `k__` will be parsed. But the rank code is not mandatory. Unassigned taxon (e.g., `s__`) and non-unique taxon are acceptable (e.g., `p__Actinobacteria` and `c__Actinobacteria`).

   Compatible with widely-used taxonomy systems in e.g., QIIME, SHOGUN, MetaPhlAn2, GTDB, etc.

5. `--map` or `-m`: Simple map of lower taxon \<tab\> higher taxon.

   One can supply multiple maps (by entering multiple `--map` parameters) to constitute several hierarchies. For example, the 1st file maps genes to UniRef entries, the 2nd maps UniRef entries to GO terms, the 3rd maps GO terms to GO slim terms, so on so forth.

   Flag `--map-as-rank` is to instruct the program to treat the map filename as rank. For example, with this flag, taxa in the 2nd column of `uniref.map.gz` will be given the rank "uniref".

Compressed files are supported and automatically recognized. For example:

```bash
woltka classify  --lineage gg_13_5_taxonomy.txt.gz ...
```


## How Woltka handles hierarchies

Hierarchy files are **additive**, i.e., if multiple files of the same or different formats are provided, all of them will be parsed and added to the classification hierarchies -- unless they conflict -- which will be noted by Woltka. Example command (example files provided under [`taxonomy`](woltka/tests/data/taxonomy), same below):

```bash
woltka classify \
  --map nucl2g.txt \
  --map g2taxid.txt \
  --nodes taxdump/nodes.dmp \
  ...
```

In this command, three layers of hierarchies are provided: 1) nucleotide ID to genome ID (`nucl2g.txt`), 2) genome ID to taxonomy ID (`g2taxid.txt`), 3) NCBI taxonomy tree (`nodes.dmp`).

**Subjects** themselves are part of the classification system. A map of subjects to one-level-higher features (e.g., a nucleotide accession to NCBI TaxID map) can be supplied with the `--map` parameter if necessary.


## Feature name dictionary

Optionally, one can supply Woltka with a feature name dictionary:

* `--names` or `-n`: NCBI-style `names.dmp` or a simple map of taxon \<tab\> name.

With a name dictionary, Woltka will append names to the output profile as an extra column. Alternatively, one may add the `--name-as-id` flag to the command and Woltka will replace feature IDs with names in the profile and the read-to-feature maps. See the [Output files](output.md) page for details.

Example command:

```bash
woltka classify \
  --input diamond.m8.gz \
  --map prot.accession2taxid.gz \
  --nodes taxdump/nodes.dmp \
  --names taxdump/names.dmp \
  --name-as-id \
  --rank genus \
  --output diamond.genus.tsv
```

Similarily, multiple name files are acceptable, and names dictionaries will be merged, unless Woltka detects any conflict among them.
