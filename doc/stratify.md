# Stratification

It is a frequent need to combine two or more classification systems in one analysis. For example, in a typical shotgun metagenomics study, one not only needs the taxonomic composition and functional potential of whole samples, but also wants to know which functional genes are found in which organism groups. This is because in many scenarios, only proteins wrapped within the same cell constitute metabolic pathways, and this information is only available when the _structure_ and _function_ of each sample are informed simultaneously.

Woltka enables this analysis through a "**stratification**" function. A stratum (plural: strata) is a subset of data grouped under the same label. Woltka allows the user to label sequences using one classification system (e.g., taxonomy), then classify sequences under each label using another classification system (e.g., function). This design maximizes flexibility and modularity. It is not bond by particular classification levels or systems, but can also work with ecological niches, phenotypic features, and other classification methods that best address specific research goals.

There are **three approaches** to achieve combined structural / functional classification. They are explained and exemplified below.


## Method I: Genome-based analysis

With this method, the only input alignments are the alignments of reads againt the whole reference **genomes**. Taxonomic classification is achieve based on the taxonomy of genomes. Meanwhile, [the alignments are translated into mapping of reads to genes based on their coordinates on the host genomes](ordinal.md). Functional classification is then made possible based on the annotation of genes.

We **recommend** this method over II (see below), for its inherent accuracy in taxonomic assignment AND observed high accuracy in functional annotation. We do note that it is computationally more expensive than the later (because there is a [translation](ordinal.md) step), but this expense is usually [affordable](perform.md).

Here is an example workflow, based on the commands and [sample data](../woltka/tests/data) introduced in the [quick-start guide](../README.md#example-usage), with merely two small modifications.

### Step 1: Taxonomic classification

```bash
woltka classify \
  -i align/bowtie2 \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank phylum \
  --name-as-id \
  --outmap mapdir \
  -o taxonomy.biom
```

The parameter `--outmap mapdir` will instruct Woltka to output read-to-taxon maps to the directory `mapdir`. These maps will serve as stratum labels for the second run:

### Step 2: Combined taxonomic/functional classification

```bash
woltka classify \
  -i align/bowtie2 \
  --coords function/coords.txt.xz \
  --map function/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --rank process \
  --stratify mapdir \
  -o taxfunc.biom
```

The parameter `--stratify mapdir` will instruct Woltka to take the previously generated taxon maps under `mapdir` to label sequences during stratification.

The output feature table, `taxfunc.biom`, is formatted like:

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


## Method II: Gene-based analysis

With this method, we start with the alignments of reads against reference **genes** (not genomes). This methods is faster, with a cost in accuracy (since genes only represent a proportion of their host genomes, and are confounded by homology, copy number, boundaries etc.). Therefore, this methods only demonstrates what "_you could_". However, if you already obtained and decided to proceed with gene alignments (for example, you already ran BLAST on UniProt), this method is worth consideration.

The following example is also based on the [sample data](../woltka/tests/data). Here we used the gene alignment results computed by DIAMOND. The gene IDs are Prodigal-annotated ORFs from the reference genomes, in the format of `genome ID <underscore> index` (e.g., `G000123456_20`). With flag `--trim-sub`, Woltka converts them to genome IDs, and everything after is straightforward.

```bash
woltka classify \
  -i align/diamond \
  --trim-sub \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank genus \
  --name-as-id \
  --outmap mapdir \
  -o taxonomy.biom
```

Then, there is no need to provide the coordinates file during functional classification:

```bash
woltka classify \
  -i align/diamond \
  --map function/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --rank process \
  --stratify mapdir \
  -o taxfunc.biom
```


## Method III: Serial collapsing

With this method, you only need to run the time-comsuming `woltka classify` workflow once (instead of twice, as in I and II) to generate a gene (ORF) profile, then you can farewell the bulky alignment files, and only work on the profile in your laptop. All you need to do then is to run as many `woltka collapse` commands as you like to collapse it to the desired structural and functional levels. This grants you the maximum flexibility. It is particularly useful when some SOP (like the [WoL SOP](wol.md) or the [Qiita server](qiita.md)) already generates the ORF profile. The results should be _identical_ to that of method I.

Step 1: Generate a gene (ORF) profile using the [coord-match algorithm](ordinal.md) without any mapping:

```bash
woltka classify -i align/bowtie2 -coords function/coords.txt.xz -o orf.biom
```

- In the resulting profile, feature IDs are like: "G000006845_1051".

Step 2 to _n_: Collapse the genome (OGU) to the desired taxonomic level:

```bash
woltka collapse -i orf.biom      -e -f 1 -m  genome-to-species.map -o species_orf.biom
woltka collapse -i species_orf.biom -f 1 -m species-to-genus.map   -o   genus_orf.biom
woltka collapse -i   genus_orf.biom -f 1 -m   genus-to-family.map  -o  family_orf.biom
...
```

- Note the first command. The `-e` or `--nested` flag breaks the ORF IDs (like "G000006845_1051") into genomes and genes, and only collapses the genome. You only need it _once_, and on the raw ORF table. Then you can stack up multiple `collapse` commands without `-e` to move up in the taxonomic hierarchy.

Step 3 to _m_: Collapse the gene (ORF) to the desired functional level:

```bash
woltka collapse -i genus_orf.biom    -f 2 -m    orf-to-ko.map      -o genus_ko.biom
woltka collapse -i genus_ko.biom     -f 2 -m     ko-to-module.map  -o genus_module.biom
woltka collapse -i genus_module.biom -f 2 -m module-to-pathway.map -o genus_pathway.biom
...
```

- This example uses the ORF by genus profile, but you can start with any intermediate profile generated from step 2. You can also mix `collapse` commands in step 2 (`-f 1`) and step 3 (`-f 2`) (order doesn't matter) until achieving the desired levels.
