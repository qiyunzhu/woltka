# Stratification

It is a frequent need to combine two or more classification systems in one analysis. For example, in a typical shotgun metagenomics study, one not only needs the taxonomic composition and functional potential of whole samples, but also wants to know which functional genes are found in which taxonomic units. This is because in many scenarios, only proteins wrapped within the same cell can constitute metabolic pathways, and this information is only available when the _taxonomy_ and _function_ of each sequence are informed simultaneously.

Woltka enables this analysis through a "**stratification**" function. A stratum (plural: strata) is a subset of data grouped under the same label. Woltka allows the user to label sequences using one classification system (e.g., taxonomy), then classify sequences under each label using another classification system (e.g., function). This design maximizes flexibility and modularity. It is not bond by particular classification levels or systems, but can also work with ecological niches, phenotypic features, and other classification methods that best address specific research goals.

Here is an example workflow, based on the commands and [sample data](../woltka/tests/data) introduced in the [quick-start guide](../README.md#example-usage), with merely two small modifications.

## Step 1: Taxonomic classification.
```
woltka classify \
  -i align/bowtie2 \
  --map taxonomy/g2tid.txt \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank phylum \
  --name-as-id \
  --outmap mapdir \
  -o taxonomy.biom
```

The parameter `--outmap mapdir` will instruct Woltka to output read-to-taxon maps to the directory `mapdir`. These maps will serve as stratum labels for the second run:

Step 2: Combined taxonomic/functional classification
```
woltka classify \
  -i align/bowtie2 \
  --coords function/coords.txt.xz \
  --map function/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --map-as-rank \
  --rank process \
  --stratify mapdir \
  -o taxfunc.biom
```

The parameter `--stratify mapdir` will instruct Woltka to take the previously generated taxon maps under `mapdir` to label sequences during stratification.

The output feature table is like:

Feature ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 |
--- | --- | --- | --- | --- |
Actinobacteria\|GO:0000003 | 10 | 6 | 2 | 0
Actinobacteria\|GO:0000005 | 4 | 20 | 3 | 0
Actinobacteria\|GO:0000006 | 0 | 12 | 5 | 2
Bacteroidetes\|GO:0000003 | 105 | 0 | 0 | 0
Bacteroidetes\|GO:0000005 | 75 | 3 | 0 | 5
Bacteroidetes\|GO:0000006 | 80 | 2 | 0 | 2
Firmicutes\|GO:0000003 | 8 | 0 | 0 | 35
Firmicutes\|GO:0000005 | 0 | 0 | 2 | 32
Firmicutes\|GO:0000006 | 0 | 0 | 0 | 18
... |
