# Tree-structured classification system

Woltka features a highly flexible hierarchical classification system. It is represented by a **tree** structure, instead of a fixed number of levels (e.g., the eight standard [taxonomic ranks](https://en.wikipedia.org/wiki/Taxonomic_rank)). It can be the traditional taxonomy, or an actual phylogenetic tree, or the hierarchies of gene-module-pathway, or any hierarchical representations of feature relationships.

The term "**rank**" (or "level") is still relevant, but it is merely a property of a feature, and does not bear information of hierarchy. For example, above a _genus_-level unit it does not have to be a _family_-level one, but could directly go to _order_, or have a _tribe_ which isn't common for the rest of the tree, or one or more nodes which do NOT have the rank assignment.

In another word, Woltka classification is **rank-independent**. This design enables finer-grain resolution of feature relationships, in addition to flexibility. It is therefore suitable for complex systems, such as phylogenetic trees.

That being said, Woltka still supports ranked hierarchies and one can instruct the program to target one or more specific ranks.


## Contents

- [Supported hierarchy files](#supported-hierarchy-files)
- [Combining hierarchies](#combining-hierarchies)
- [Taxon name dictionary](#taxon-name-dictionary)
- [How Woltka handles hierarchies](#how-woltka-handles-hierarchies)
- [Hierarchy file specs](#hierarchy-file-specifications)
- [Multiple mapping](#multiple-mapping)


## Supported hierarchy files

Wolkta supports various types and formats of classification systems, as listed below. Details specs of formats are provided at the [bottom](#hierarchy-file-specifications) of this page.

1. `--nodes`: [NCBI-style](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) `nodes.dmp` (columns are delimited by "\<tab\>|\<tab\>") or a tab-delimited file, in which each taxon (1st column) points to its parent taxon (2nd column). Rank (3rd column) is optional.

2. `--newick`: [Newick-format](https://en.wikipedia.org/wiki/Newick_format) tree, in which labels of nodes (tips, internal nodes and root) are considered as taxa. All nodes must have labels and all labels must be unique.

3. `--lineage`: Map of taxon to lineage string (`;`-delimited taxa from high to low)

   It can be [Greengenes-style](http://greengenes.secondgenome.com/) taxonomy where rank codes such as `k__` will be parsed. But the rank code is not mandatory. Unassigned taxon (e.g., `s__`) and non-unique taxon are acceptable (e.g., `p__Actinobacteria` and `c__Actinobacteria`).

   Compatible with widely-used taxonomy systems in e.g., QIIME, SHOGUN, MetaPhlAn2, GTDB, etc.

3. `--columns`: Table of per-rank assignments. Each column represents a rank. The column header will be treated as the rank name.

5. `--map` or `-m`: Simple map of lower taxon \<tab\> higher taxon.

   Flag `--map-as-rank` is to instruct the program to treat the map filename as rank. For example, with this flag, taxa in the 2nd column of `uniref.map.gz` will be given the rank "uniref".

Compressed files are supported and automatically recognized. For example, reading the gzipped Greengenes taxonomy file is as simple as:

```bash
woltka classify --lineage gg_13_5_taxonomy.txt.gz ...
```

## Combining hierarchies

Hierarchy files are **additive**, i.e., if multiple files of the same or different formats are provided, all of them will be parsed and added to the classification hierarchies -- unless they conflict -- which will be noted by Woltka. Example command (example files provided under [`taxonomy`](woltka/tests/data/taxonomy), same below):

```bash
woltka classify \
  --map nucl2g.txt \
  --map g2taxid.txt \
  --nodes nodes.dmp \
  ...
```

In this command, three layers of hierarchies are provided: 1) nucleotide ID to genome ID (`nucl2g.txt`), 2) genome ID to taxonomy ID (`g2taxid.txt`), 3) NCBI taxonomy tree (`nodes.dmp`).

Another example, in which three simple mapping files are provided, allowing Woltka to group genes into UniRef entries, then to GO terms, then to GO slim terms.

```bash
woltka classify \
  --map gene2uniref.map \
  --map uniref2go.map \
  --map go2goslim.map \
  --rank go,goslim \
  ...
```

A third example, in which two [GTDB](gtdb.md) lineage files (Bacteria and Archaea) are directly read without the need for merging them beforehand:

```bash
woltka classify \
  --lineage bac120_taxonomy.tsv \
  --lineage ar122_taxonomy.tsv \
  --rank species \
  ...
```


## Taxon name dictionary

Optionally, one can supply Woltka with a taxon (a.k.a. classification unit, or feature) name dictionary using the `--names` (or `-n`) parameter. This file can be an NCBI-style `names.dmp` (example: [taxonomy/names.dmp](../woltka/tests/data/taxonomy/names.dmp) or a simple map of `taxon <tab> name`.

With a name dictionary, Woltka will append names to the output profile as an extra column. Alternatively, one may add the `--name-as-id` flag to the command and Woltka will replace feature (taxon) IDs with names in the profile and the read-to-feature maps. See the [Output files](output.md) page for details.

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


## How Woltka handles hierarchies

The data structure of the classification hierarchies in Woltka is a mapping of taxon to its direct parent. It is a simple **tree** structure, allowing the program to move from any taxon up in the hierarchy until it reaches the root.

The tree always has a single **root**. When the supplied classification file does not have a unique top-level taxon, or multiple classifcation files are supplied which do not have overlaps, Woltka will automatically create a root and assign all top-level taxa as its child nodes.

Furthermore, **subjects** themselves are part of the classification system. Therefore in Woltka, the classification hierarchy is a massive tree from all subjects to a single root.


## Hierarchy file specifications

### 1. Nodes (`--nodes`)

The classical format used by the NCBI taxonomy database (i.e., "[taxdump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)"). It is a multi-column text file with columns delimited by "\<tab\>|\<tab\>".

Woltka cares only the first three columns: 1st: the **taxon** (i.e., TaxID), 2nd: its **parent** taxon, and 3rd: the **rank**. The 3rd column is optional, because Woltka uses but [does not rely](classify.md) on ranks.

Therefore, it is okay to provide a **custom** nodes file, with only first three or two columns.

To make things simple, Woltka further allows using a simple \<tab\> as the column delimiter.

Example (simplified from [taxonomy/nodes.dmp](../woltka/tests/data/taxonomy/nodes.dmp)):

```
1 <tab> 1 <tab> root
2 <tab> 1 <tab> domain
2157 <tab> 1 <tab> domain
1224 <tab> 2 <tab> phylum
...
```

### 2. Newick tree (`--newick`)

[Newick](https://en.wikipedia.org/wiki/Newick_format) is the standard file format for phylogenetic trees. Woltka supports direct conversion of Newick trees into classification hierarchies, in which labels of nodes (tips, internal nodes and root) are considered as taxa. _All nodes must have labels and all labels must be unique_.

For example, a simple Newick string `((a,b),(c,d));` is not valid for Woltka, because it lacks internal node labels. It needs to be `((a,b)c,(d,e)f)g;` instead.

In this simple example, we obtained seven taxa. We see that `c` is the parent of `a` and `b`, and `g` is the parent of `c` and `f`, etc.

Branch lengths (e.g., `((a:1.2,b:0.8)c:0.3)`) are tolerated but will be stripped in Woltka. That is, Woltka only considers **topology** but not phylogenetic distances.

In the Web of Life ([WoL](https://biocore.github.io/wol/)) database, the [reference phylogenetic tree](https://biocore.github.io/wol/data/trees/astral/branch_length/cons/collapsed/astral.cons.nid.e5p50.nwk) already has unique internal node labels (`N###`) in addition to tip labels (genome IDs, `G###`).

Example (refer to [tree.nwk](../woltka/tests/data/tree.nwk)):

```
((((((G000191145,G000012245)N63,((((G000195995,G000008865)N94,(G000025565,G000240185)N95)N91,G000009065)N84,...
```

### 3. Lineage strings (`--lineage`)

A **lineage string** is the concatenation of all taxa (classification units) from most ancestral to current taxon, delimited by a semicolon (`;`), with or without a trailing space.

```
Organic compounds; Lipids; Steroids; Bile acids; Glycinated bile acids
```

Woltka reads a mapping of `subject <tab> lineage` and constructs hierarchies.

In taxonomic classification, a common practice is to add a single-letter prefix to each taxon, representing its taxonomic rank. For example (see [taxonomy/lineage.txt](../woltka/tests/data/taxonomy/lineage.txt)):

```
Seq1 <tab> k__Bacteria; p__Firmicutes; c__Bacilli
Seq2 <tab> k__Bacteria; p__Firmicutes; c__Clostridia
Seq3 <tab> k__Bacteria; p__Tenericutes; c__Mollicutes
```

This format, a.k.a., [Greengenes-style](http://greengenes.secondgenome.com/) lineage, is widely used in bioinformatics tools and databases, such as QIIME, SHOGUN, MetaPhlAn2 and GTDB.

Woltka automatically recognizes the rank codes and convert them to rank names. The eight standard ranks are supported: **k**ingdom; **p**hylum; **c**lass; **o**rder, **f**amily, **g**enus, **s**pecies, and s**t**rain.

In addition, (note!) **d**omain will be parsed as kingdom. You may disagree with this confusion (as I do), but it is an operational solution compatible with existing protocols.

Unlike other types of hierarchies, with lineages, Woltka reports the _entire lineage string_ in the output. For example, instead of `c__Bacilli`, Woltka will report `k__Bacteria; p__Firmicutes; c__Bacilli`. This is for compatibility with existing tools (such as QIIME 2). So there is no need to explicitly append the lineage string (`--add-lineage`).

Empty levels in the end are discarded. e.g., `k__Bacteria; p__` is not a valid taxon. Empty levels in the middle are kept. e.g., `k__Bacteria; p__; c_Clostridia` will not be shortened into `k__Bacteria;c_Clostridia`.

[Note] With many other bioinformatics tools, the use of lineage implicates a fixed number of ranks for all taxa. However in Woltka, lineages are still treated in a rank-flexible manner. For example, the following [lineage](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=9606&lvl=3&lin=f&keep=1&srchmode=1&unlock) is valid:

```
cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; Homininae; Homo; Homo sapiens
```

### 4. Table of per-rank columns (`--columns`)

A tab-delimited table (TSV), with aach column representing a rank, and the column header serving as the rank name. For example ([taxonomy/rank_names.tsv](../woltka/tests/data/taxonomy/rank_names.tsv)):

```
ID <tab> kingdom <tab> phylum <tab> class ...
Seq1 <tab> Bacteria <tab> Proteobacteria <tab> Gammaproteobacteria
Seq2 <tab> Bacteria <tab> Firmicutes <tab> Bacilli
Seq3 <tab> Bacteria <tab> Proteobacteria <tab> Betaproteobacteria
...
```

The assumption with such as table is that ranks displayed in the header are always from high to low. Therefore, it is a form of **fixed-rank** classification.

Woltka will construct classification hierarchies from this table, and will report if it encounters any conflicting relationships among levels.

### 5. Simple map (`--map` or `-m`)

The simplest form of mapping of `taxon <tab> parent`. For example ([taxonomy/g2tid.txt](../woltka/tests/data/taxonomy/g2tid.txt)):

```
G000006745 <tab> 243277
G000006785 <tab> 160490
G000006845 <tab> 242231
...
```

Only the first two columns are considered, if there are more of them.

With flag `--map-as-rank`, Woltka will extract a **rank** name from the filename to represent the **parents** (2nd column). The naming rule is like follows:

- `uniref.map` => `uniref`
- `prot2taxid.txt.gz` => `taxid`
- `reaction_to_pathway.tsv` => `pathway`
- `apple-to-orange` => `orange`


## Multiple mapping

As explained above, Woltka's main classification workflow relies on a tree-structured hierarchy, i.e., a classification unit can only point to one parent classification unit. However, in some scenarios one may want to perform one-to-many mapping. Examples include functional analyses, where one gene may be involved in multiple pathways. Woltka provides a [profile collapsing](collapse.md) tool to enable this analysis ([see details](collapse.md)).
