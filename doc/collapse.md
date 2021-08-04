# Profile collapsing

The **profile collapsing** function is a lightweight and flexible addition to the main classification workflow. It allows the user to convert an _existing_ profile based on a mapping of source features to target features. It highlights the support for **many-to-many mapping**.

```bash
woltka tools collapse -i input.biom -m mapping.txt -o output.biom
```


## Contents

- [Use cases](#use-cases)
- [Mapping file format](#mapping-file-format)
- [Parameters](#parameters)
- [Considerations](#considerations)


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


## Mapping file format

A mapping file is a text file with entries separated by tabs. The number of fields per line is _arbitrary_. The first field is the source feature ID. The second to last fields are target feature ID(s). Duplicates in sources or targets are allowed. For examples:

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

## Parameters

### Division

By default, if one source feature is simultaneously mapped to _k_ targets, each target will be counted once. With the `--divide` or `-d` flag added to the command, each target will be counted 1 / _k_ times.

Whether to enable division depends on the nature and aim of the analysis. For example, one gene is involved in two metabolic pathways (which isn't uncommon), should each pathway be counted once, or half time?

### Stratification

Woltka supports collapsing a [stratified](stratify.md) profile using one field in the feature IDs. This can be done using the `--field` or `-f` parameter followed by the field index (starting from 1). For example, if the feature IDs are in the format of "Species|Gene", one may collapse genes into pathways while keeping species by adding `-f 2`.

### Feature names

Once a profile is collapsed, the metadata of the source features ("Name", "Rank", and "Lineage") will be discarded. One may choose to supply a target feature name file by `--names` or `-n`, which will instruct the program to append names to the profile as a metadata column ("Name").


## Considerations

It is important to highlight that one-to-many mapping may change some of the
underlying statistical assumptions of downstream analyses.

In the default mode, because one source may be collapsed into multiple targets, the total feature count per sample may be inflated, and the relative abundance of each feature may no longer correspond to that of the sequences assigned to it. In other words, this breaks the [compositionality](https://en.wikipedia.org/wiki/Compositional_data) of the data.

How significantly this may impact an analysis depends on the relative frequency of multiple mappings found in the data, the biological relevance of the affected features, and the statistical nature of the analysis.

For example, in the `reaction-to-ec.txt` file under [MetaCyc](metacyc.md), 80 out of 3618 (2.2%) reactions have more than one corresponding EC number. Whether such a translation may be considered as unique (and whether the resulting table is still compositional) is a call of the user.

A solution to this is to turn on the [division](#division) flag (`-d`). This guarantees that the sum of feature counts remains the same after collapsing. But one should consider the biological implication before making a decision (see [above](#division)).
