# Profile collapsing

The **profile collapsing** function is a lightweight and flexible addition to the main classification workflow. It allows the user to convert an _existing_ profile based on a mapping of source features to target features, and/or internal hierarchies of features. It highlights the support for **many-to-many mapping**. For examples:

Collapse features based on a mapping file:

```bash
woltka collapse -i input.biom -m mapping.txt -o output.biom
```

Collapse nested features to the first level:

```bash
woltka collapse -i input.biom -e -f 1 -o output.biom
```


## Contents

- [Use cases](#use-cases)
- [Feature mapping](#feature-mapping)
  - [Mapping file](#mapping-file)
  - [Value division](#value-division)
  - [Considerations](#considerations)
- [Stratified features](#stratified-features)
  - [Collapse to field](#collapse-to-field)
  - [Collapse by mapping](#collapse-field-by-mapping)
- [Nested features](#nested-features)
  - [Collapse to level](#collapse-to-level)
  - [Collapse by mapping](#collapse-level-by-mapping)


## Use cases

With this tool one can achieve the following goals:

1. Translate feature IDs into names or descriptions. Examples:
   - Translate taxonomic IDs to taxon names.
   - Translate ORF IDs to gene IDs, while **dropping** the unannotated.
   - Translate [UniRef](https://www.uniprot.org/help/uniref) IDs to protein names, while **merging** same names.

2. Group lower features into higher categories. Examples:
   - Group genera into families, then into orders.
   - Group chemical structures by chemical ontology.

3. Convert lower features into higher ones, where each lower feature may correspond to **multiple** higher features. Examples:
   - Convert KEGG [orthologs](https://www.genome.jp/kegg/ko.html) to [pathways](https://www.genome.jp/kegg/pathway.html).
   - Convert [GO](http://geneontology.org/docs/ontology-documentation/) terms to [GO Slim](http://www-legacy.geneontology.org/GO.slims.shtml) terms.

The last usage is an important complement to the main classification workflow, which currently relies on a tree structure and does not support one-to-many mapping. This can be achieved by using the profile collapsing function (although one can only move up one level per run).

See [considerations](#considerations) below for a discussion of the potential change of statistical properties of data.


## Feature mapping

The basic usage of the `collapse` command is to translate feature IDs using an external mapping file.


### Mapping file

A mapping file (specified by `--map` or `-m`) is a text file with entries separated by tabs. The number of fields per line is _arbitrary_. The first field is the source feature ID. The second to last fields are target feature ID(s). Duplicates in sources or targets are allowed. For examples:

1. One/many-to-one:
```
source1 <tab> target1
source2 <tab> target2
source3 <tab> target2
source4 <tab> target3
...
```

2. Many-to-many (multiple targets per line):
```
source1 <tab> target1
source2 <tab> target1 <tab> target2
source3 <tab> target2 <tab> target3 <tab> target4
...
```

3. Many-to-many (multiple same sources):
```
source1 <tab> target1
source1 <tab> target2
source2 <tab> target2
source3 <tab> target3
source4 <tab> target3
...
```

Once a profile is collapsed, the metadata of the source features ("Name", "Rank", and "Lineage") will be discarded. One may choose to supply a target feature name file by `--names` or `-n`, which will instruct the program to append names to the profile as a metadata column ("Name").

### Value division

By default, if one source feature is simultaneously mapped to _k_ targets, each target will be counted once. With the `--divide` or `-d` flag added to the command, each target will be counted 1 / _k_ times.

Whether to enable division depends on the nature and aim of the analysis. For example, one gene is involved in two metabolic pathways (which isn't uncommon), should each pathway be counted once, or half time? The user needs to make a decision.

### Considerations

It is important to note that one-to-many mapping may change some of the underlying statistical assumptions of downstream analyses.

In the default mode, because one source may be collapsed into multiple targets, the total feature count per sample may be inflated, and the relative abundance of each feature may no longer correspond to that of the sequences assigned to it. In other words, this breaks the [compositionality](https://en.wikipedia.org/wiki/Compositional_data) of the data.

How significantly this may impact an analysis depends on the relative frequency of multiple mappings found in the data, the biological relevance of the affected features, and the statistical nature of the analysis.

For example, in the `reaction-to-ec.txt` file under [MetaCyc](metacyc.md), 80 out of 3618 (2.2%) reactions have more than one corresponding EC number. Whether such a translation may be considered as unique (and whether the resulting table is still compositional) is a call of the user.

A solution to this is to turn on the [division](#division) flag (`-d`). This guarantees that the sum of feature counts remains the same after collapsing. But one should consider the biological implication before making a decision (see [above](#division)).


## Stratified features

Woltka supports collapsing a [stratified](stratify.md) profile using one field in the feature IDs. This can be done using the `--field` or `-f` parameter followed by the field index (starting from 1). The default field delimiter is a pipe (`|`), but one can customize it using `--sep` or `-s`.

For example, in the following profile `species_gene.tsv`, feature IDs has the format of "species|gene", representing particular genes (e.g., _rpoA_) found in particular species (e.g., _E. coli_).

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`Ecoli\|rpoA` | 4 | 15 | 0
`Ecoli\|rpoB` | 12 | 7 | 5
`Sente\|rpoC` | 9 | 0 | 3
`Cdiff\|ftsZ` | 1 | 6 | 0

### Collapse to field

One can collapse the "species|gene" features into just species (regardless of gene) with:

```bash
woltka collapse -i species_gene.tsv -f 1 -o species.tsv
```

The output profile `species.tsv` is like:

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`Ecoli` | 16 | 22 | 5
`Sente` | 9 | 0 | 3
`Cdiff` | 1 | 6 | 0

Alternatively, one can use `-f 2` to collapse to genes (regardless of species).

### Collapse field by mapping

With a species-to-phylum mapping file `phylum.map`, one can collapse the first field (species):

```bash
woltka collapse -i species_gene.tsv -f 1 -m phylum.map -o phylum_gene.tsv
```

The output profile `phylum_gene.tsv` will be like:

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`Proteo\|rpoA` | 4 | 15 | 0
`Proteo\|rpoB` | 12 | 7 | 5
`Proteo\|rpoC` | 9 | 0 | 3
`Firmic\|ftsZ` | 1 | 6 | 0

With a gene-to-biological process mapping file `process.map`, one can collapse the second field (gene):

```bash
woltka collapse -i species_gene.tsv -f 2 -m process.map -o species_process.tsv
```

The output profile `species_process.tsv` will be like:

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`Ecoli\|mRNASyn` | 16 | 22 | 5
`Sente\|mRNASyn` | 9 | 0 | 3
`Cdiff\|CellDiv` | 1 | 6 | 0

One can combine the two operations:

```bash
woltka collapse -i species_gene.tsv -f 1 -m phylum.map -o phylum_gene.tsv
woltka collapse -i phylum_gene.tsv -f 2 -m process.map -o phylum_process.tsv
```

The output profile `phylum_process.tsv` will be like:

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`Proteo\|mRNASyn` | 25 | 22 | 8
`Firmic\|CellDiv` | 1 | 6 | 0


## Nested features

In some scenarios, a feature ID itself contains hierarchical information. This is similar to a stratified feature (see above) in that it contains multiple fields, but each field is a child of the previous field, and it only makes sense when written after the previous one. We refer to them as **nested features**.

For example, "G12_34" represents the 34th ORF annotated from genome 12. In other words, this is equivalent to a stratified feature ID "G12_34|G12_34", in which the 1st field represents the genome ([OGU](ogu.md)) and the 2nd represents the gene (ORF). Such profiles can be generated using Woltka's ["coord-match" functional profiling](ordinal.md) function.

The `--nested` or `-e` flag instructs Woltka to treat feature IDs as nested. The default delimiter for identifying fields is an underscore (`_`) (in contrast to a pipe (`|`) in a stratified feature), but one may customize it using `-s` (see above). Then, with the `--field` or `-f` parameter (see above), one can specify the level in the nested features on which collapsing will occur.

### Collapse to level

For example, with an ORF table `orf.tsv`, one can do:

```bash
woltka collapse -i orf.tsv -e -f 1 -o ogu.tsv
```

This will collapse an ORF table into an [OGU table](ogu.md), in which only "G12" but not "_34" is retained.

- Note howerver, that the resulting OGU table may not be identical to the one produced using the [OGU protocol](ogu.md), which also considers intergenic regions (not just coding genes) in a genome.

In this example, one can also do `-f 2`, but nothing will change (because there are two levels in total.)

In the following example `ec4.tsv`, feature IDs are 4-level [EC numbers](https://en.wikipedia.org/wiki/Enzyme_Commission_number):

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`EC:6.4.1.2` | 36 | 8 | 22
`EC:6.4.1.3` | 9 | 13 | 14
`EC:6.3.4.14` | 62 | 4 | 20
`EC:2.3.1.85` | 5 | 16 | 11

One can collapse them into 2-level EC numbers with:

```bash
woltka collapse -i ec4.tsv -e -s "." -f 2 -o ec2.tsv
```

The output profile `ec2.tsv` is like:

Feature ID | Sample 1 | Sample 2 | Sample 3
--- | --- | --- | ---
`EC:6.4` | 45 | 21 | 36
`EC:6.3` | 62 | 4 | 20
`EC:2.3` | 5 | 16 | 11


The following example `free.tsv` is a taxonomic profile with arbitrary ranks. It can be generated using Woltka's [free-rank classification](classify.md#target-rank-or-no-rank), or other programs such as QIIME 2 and MetaPhlAn.

Feature ID | Sample 1 | Sample 2
--- | --- | ---
`p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus` | 2 | 0
`p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia` | 8 | 3
`p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae` | 37 | 16
`p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__fragilis` | 56 | 24
`p__Proteobacteria; c__Epsilonproteobacteria; o__Campylobacterales` | 0 | 7
`p__Spirochaetes` | 13 | 9

One can collapse them into the class level with:

```bash
woltka collapse -i free.tsv -e -s "; " -f 2 -o class.tsv
```

The output profile `class.tsv` is like:

Feature ID | Sample 1 | Sample 2
--- | --- | ---
`p__Firmicutes; c__Bacilli` | 2 | 0
`p__Firmicutes; c__Clostridia` | 45 | 19
`p__Bacteroidetes; c__Bacteroidia` | 56 | 24
`p__Proteobacteria; c__Epsilonproteobacteria` | 0 | 7

Note that nested features do not have to have the same number of levels. Those with fewer levels than the specified field index will be dropped.

### Collapse level by mapping

One can collapse a certain level in nested features (`-e`) using a mapping file (`-m`).

For example, in an ORF table `orf.tsv` (see [details](ordinal.md)), feature IDs are like:

```
G000007845_1274
```

With a genome-to-genus mapping file `genus.map`, one can collapse the first field (i.e., taxonomic analysis):

```bash
woltka collapse -i orf.biom -e -f 1 -m genus.map -o genus_orf.biom
```

The resulting feature IDs are like:

```
Bacillus|G000006785_1274
```

With an ORF-to-UniRef mapping file `uniref.map`, one can collapse the second field (i.e., functional analysis):

```bash
woltka collapse -i orf.biom -e -f 2 -m uniref.map -o ogu_uniref.biom
```

The resulting feature IDs are like:

```
G000006785|J9GI19
```

One can execute the `collapse` command multiple times to collapse the first and second fields separately to achieve desired taxonomic and functional resolution. Note that starting from the second command, features are no longer nested.

```bash
woltka collapse -i orf.biom -e -f 2 -m uniref.map -o ogu_uniref.biom
woltka collapse -i ogu_uniref.biom -f 1 -m genus.map -o genus_uniref.biom
woltka collapse -i genus_uniref.biom -f 2 -m uniref2go.map -o genus_go.biom
...
```

Note that feature IDs in the original profile are nested. However, in the collapsed profile, they become stratified by a pipe (`|`) from the collapsed field to the end. This is because the collapsed field has broken the nestedness.

It is important to avoid confusion when collapsing a certain level in nested features. For example, the original feature IDs are like (4-level EC numbers):

```
EC:2.7.2.10
```

One can collapse the 2nd level by using an enzyme-to-substrate mapping file `substrate.map`:

```bash
woltka collapse -i ec4.tsv -e -f 2 -m substrate.map -o output.tsv
```

The output feature IDs will be like:

```
EC:2|phosphate|EC:2.7.2.10
```
