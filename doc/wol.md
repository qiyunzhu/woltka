# Working with "Web of Life" (WoL)

**TL;DR**: Here is a script for running the full-board WoL/Woltka workflow, including genome (OGU) and gene profilings, taxonomic classification (using NCBI taxonomy), and functional classification (using UniRef, GO, MetaCyc and optionally KEGG). You need to download the WoL data release, perform sequence alignment, customize the "Parameters" section of the script, then run it.

- [**wolsop.sh**](wolsop.sh) (click to download)

Meanwhile, a basic WoL/Woltka workflow is available from our web-based microbiome study platform: **Qiita** (https://qiita.ucsd.edu/) (see [details](qiita.md)).

## The WoL database

The "Web of Life" (WoL) project aims to reconstruct an accurate reference phylogeny for microbial genomes, and to build resources that can benefit microbiome researchers. In phase I ([Zhu et al., 2019](https://www.nature.com/articles/s41467-019-13443-4)), we built a reference [tree](https://biocore.github.io/wol/data/trees/tree.nwk) of 10,575 bacterial and archaeal genomes using 381 marker genes. The basic WoL data release, which contains all necessary files for performing microbiome data analysis, is available for download from this FTP site:

```bash
wget -r -np -nH --cut-dirs=2 ftp://ftp.microbio.me/pub/wol-20April2021/
```

Alternatively, you may download the database from the following Dropbox site. Specifically, `wolr1_seqs` is the genome sequences with which you can build a database for sequence alignment, which will happen prior to running Woltka; `wolr1_maps` contains all mapping files which you will need for doing taxonomic, functional and phylogenetic analyses using Woltka.

- https://www.dropbox.com/sh/rwhetdmwomi9k34/AAAKRLIJ2nFniJESgb8lHR-na?dl=0

Meanwhile, the full WoL data release, which also contains raw and processed sequence data, metadata, alternative trees, and pre-built databases for multiple metagenomics tools, are available from our Globus endpoint: [WebOfLife](https://app.globus.org/file-manager/collections/31acbeb8-c62f-11ea-bef9-0e716405a293) (see [instruction](https://biocore.github.io/wol/download#download-via-globus)). The project is detailed at our website: https://biocore.github.io/wol/, including documentation, code, protocols, a gallery and a visualizer.

The following tutorial assumes that you have downloaded the basic or full WoL data release directory. The paths mentioned below are relative to this directory. Specifically, the following directories and files are relevant:

- `databases/bowtie2/`
- `proteins/coords.txt.xz`
- `taxonomy/`
- `function/`


## Database indexing

The full WoL data release contains a pre-built [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) index which can be directly used (and you can skip this section). The basic release only contains the original genome sequences (`concat.fna.xz`). You will need to build a Bowtie2 index by yourself (see [here](align.md#database-indexing) for details):

```bash
mkdir -p databases/bowtie2
xz -d -k concat.fna.xz
bowtie2-build --threads 8 concat.fna databases/bowtie2/WoLr1
rm concat.fna
```

## Sequence alignment

With an indexed database, you can align your sequencing data (e.g., paired-end sequences `R1.fq` and `R2.fq`) against reference genomes using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (see [here](align.md#sequence-alignment) for details):

```bash
bowtie2 -p 8 -x db -1 R1.fq -2 R2.fq -S output.sam --very-sensitive --no-head --no-unal
```

This step will generate a [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) format alignment file.

**Note**: You can also run Bowtie2 manually using your choice of parameters, or using other aligners and corresponding databases. Woltka is designed for flexibility.


## The OGU analysis

**OGU** (operational genomic unit) ([Zhu et al., 2022](https://journals.asm.org/doi/10.1128/msystems.00167-22)) is a notion we proposed to define the minimum unit of microbiome composition allowed by shotgun metagenomic data. OGUs are reference genomes to which any input sequences have matches. This maximizes the resolution of microbiome composition, and allows for phylogeny-aware analyses using the WoL reference phylogeny. See [here](ogu.md) for details.

In Woltka, an OGU analysis is simply a classification process without a classification system. In such case, the output features will just be the subjects in the alignment file, namely reference genomes.

```bash
woltka classify -i input.sam -o ogu.biom
```

Woltka accepts a multiplexed alignment file, or a directory of multiple per-sample files, or a mapping of sample IDs to filepaths as input. See [here](input.md#input-filepath) for details. Woltka automatically parses compressed alignment files (such as `input.sam.gz`) so you may [compress](perform.md#compress-alignment-files) them to save disk space.

The resulting file, `ogu.biom` is a **feature table** describing the composition of the microbiome using OGUs (reference genome IDs, such as `G000123456`). This file can be analyzed using standard microbiome data analysis packages such as [QIIME 2](https://qiime2.org/), in combination with the WoL tree. Here is an [example](ogu.md#ogu-analysis-using-qiime-2).


## Free-rank classification

With a tree-structured classification system (e.g., taxonomy), Woltka will assign each query sequence to the lowest classification unit that can describe all its subjects, which can be one exact match, or multiple matches that are equally or similarily well aligned. The resulting features are not restricted to a specific taxonomic rank. We refer to it as "**free-rank**" classification.

**Note**: In the resulting feature table, features may be nested (e.g., both _Escherichia_ and _Eschierichia coli_ are present in the same table), but their **frequencies are exclusive**. Therefore, the compositionality of the table is retained. Please note this distinction from outputs of other programs such as MetaPhlAn and Kraken.

The following command performs classification on the original NCBI taxonomy (taxdump):

```bash
woltka classify \
  --input  input.sam \
  --map    taxonomy/taxid.map \
  --nodes  taxonomy/nodes.dmp \
  --names  taxonomy/names.dmp \
  --output free.biom
```

- The `--names` parameter can be omitted if you don't need names or there are no names available. Same below.

Alternatively, one may use lineage strings extracted from NCBI (will lose some resolution, but results are more structured, especially for users familiar with QIIME 2):

```bash
woltka classify \
  --input   input.sam \
  --lineage taxonomy/lineages.txt \
  --output  free.biom
```

We also provide original and curated [NCBI](https://biocore.github.io/wol/data/taxonomy/ncbi/) and [GTDB](https://biocore.github.io/wol/data/taxonomy/gtdb/) taxonomy for choice.


## Classification by taxonomic ranks

One may specific one or multiple classification ranks in the command, and Woltka will generate one feature table for each rank, in which all features are classification units at this rank. This method is inline with the traditional **taxonomic profiling**.

```bash
woltka classify \
  --input   input.sam \
  --lineage taxonomy/lineages.txt \
  --rank    phylum,genus,species \
  --output  .
```

You will get three output files: `phylum.biom`, `genus.biom` and `species.biom`.

**Note**: There three files are automatically generated in Qiita's Woltka workflow (see above).


## "Coord-match" gene profiling

Woltka implements an efficient algorithm which assigns query sequences to open reading frames (ORFs) of the reference genomes based on coordinates matching. See [details](ordinal.md). This new feature combines taxonomic and functional analyses of the microbiome in one alignment step, ensuring consistency and also improving accuracy (manuscript in preparation).

```bash
woltka classify \
  --input  input.sam \
  --coords proteins/coords.txt.xz \
  --output gene.biom
```

The file `coords.txt.xz` defines the coordinates of ORFs on reference genomes (again, compressed files are supported). In the output feature table, each feature is an ORF ID (like `G000123456_789`, meaning the 789th ORF on genome G000123456).

While informing the genetic composition of the microbiomes, this table itself may not be suitable for community analyses (because it is too sparse). But it serves as the basis for stacking up higher-level functional categories, as explained below.


## Functional classification of genes

Woltka provides flexible supports for various classification systems, including ones that inform the functional capacities of organisms. There are two methods for functional classification using Woltka.

### Method 1

Combine gene profiling and one or multiple higher-level functional profilings in one `classify` command.

```bash
woltka classify \
  --input  input.sam \
  --coords proteins/coords.txt.xz \
  --map    function/uniref/uniref.map.xz \
  --names  function/uniref/uniref.name.xz \
  --map    function/kegg/ko.map.xz \
  --names  function/kegg/ko.name \
  --rank   uniref,ko \
  --to-tsv \
  --output .
```

In this command, `uniref.map` is a mapping of ORFs to UniRef entries. `ko.map` is a mapping of UniRef entries to KEGG Ontologies (KOs). The two `.name` files provide descriptive names of the UniRef and KO IDs, which will be appended to the tables.

The output files will be `uniref.tsv` and `ko.tsv`. The optional flag `--to-tsv` instructs the program to generate tab-delimited plain text tables instead BIOM format. This applies to taxonomic profiling too. See [details](output.md).

- It is useful when you (human) want to read the names of the entries. But even if the default BIOM format is used, the names will still be appended as a metadata column (not currently supported by QIIME 2 though) See [details](output.md#biom).

### Method 2

Starting from the gene profile (`gene.biom`), sequentially stack up functional hierarchies, one at a time, using the `collapse` command (see [details](collapse.md)).

The following commands generate profiles describing the metabolic capacities of the microbiome at multiple levels defined in [MetaCyc](https://metacyc.org/) (see [details](metacyc.md)).

Make an abbreviation:

```bash
mcdir=function/metacyc
```

Collapse ORFs into MetaCyc proteins:

```bash
woltka collapse -i gene.biom -m $mc/protein.map.xz -n $mc/protein_name.txt -o protein.biom
```
- Here `-i`, `-o`, `-m` and `-n` are abbreviations for `--input`, `--output`, `--map` and `--names`.

Collapse proteins into enzymatic reactions:

```bash
woltka collapse -i protein.biom -m $mc/protein-to-enzrxn.txt -n $mc/enzrxn_name.txt -o enzrxn.biom
```

Collapse enzymatic reactions into reactions:

```bash
woltka collapse -i enzrxn.biom -m $mc/enzrxn-to-reaction.txt -n $mc/reaction_name.txt -o reaction.biom
```

Collapse reactions to metabolic pathways:

```bash
woltka collapse -i reaction.biom -m $mc/reaction-to-pathway.txt -n $mc/pathway_name.txt -o pathway.biom
```

So on so forth. See [here](metacyc.md) for a graph of all available collapsing directions.

### Comparison

**Important**: The differences between [method 1](#method-1) (`classify`) and [method 2](#method-2) (`collapse`) are:

`classify` only supports a tree structure, in which one child unit has exactly one parent unit. This is typical in taxonomic classification. If multiple parents are present, all but the first parent will be discarded. In contrast, `collapse` supports **one-to-multiple** mappings, therefore it is more suitable when this is the norm instead of exception, especially in functional classification (where one gene can be involved in multiple metabolic pathways).

`classify` always ensures the [**compositionality**](https://en.wikipedia.org/wiki/Compositional_data) of the feature table, in which the frequencies match the numbers of aligned sequences. `collapse` however does not by default: In a one-to-multiple mapping, all parents will be counted once. But one can add `-d` to the `collapse` command to divide the counts by the number of parents so that the compositionality is retained. See [here](collapse.md#considerations) for a detailed discussion.

However, if the mappings are **unique** (one-to-one), the two methods produce mutually identical results (minor differences may arise during number rounding), and the concern of compositionality is no longer relevant.

Here is a list of currently released primary and secondary mappings and their uniqueness:

Source | Target | Unique
--- | --- | ---
WoL protein | UniRef | Yes
WoL protein | MetaCyc protein | Yes
UniRef | KEGG ontology (KO) | No
UniRef | Gene Ontology (GO) | No
UniRef | OrthoDB | Yes
UniRef | EggNOG | No
UniRef | RefSeq | No

Meanwhile, all name files (such as `gene_name.txt` under MetaCyc) are unique.


## Stratified taxonomic / functional classification

Woltka allows combining taxonomic and functional classifications to generate **stratified** results. For example, you want to get the composition of KEGG ontologies of each genus:

First, run **taxonomic** classification at the genus level, and export read-to-genus maps:

```bash
woltka classify \
  --input  input.sam \
  --map    taxonomy/taxid.map \
  --nodes  taxonomy/nodes.dmp \
  --names  taxonomy/names.dmp \
  --rank   genus \
  --name-as-id \
  --output genus.biom \
  --outmap mapdir
```

This command generates mappings of sequences to genera, in addition to a genus-level feature table. This information will guide the subsequent stratification.

Second, run **functional** annotation, with the read-to-genus maps incorporated:

```bash
woltka classify \
  --input    input.sam \
  --coords   proteins/coords.txt.xz \
  --map      function/uniref/uniref.map.xz \
  --map      function/kegg/ko.map.xz \
  --rank     ko \
  --stratify mapdir \
  --output   ko_by_genus.biom
```

In the output table, features will be like `Escherichia|K00133`, `Salmonella|K00604`, etc. See [here](stratify.md) for mode details about stratification.

With a stratified taxonomic/functional profile, one may still perform further [collapsing](#method-2) using one of the two systems. For example:

```bash
woltka collapse -f 2 -i ko_by_genus.biom -m ko/ko-to-go.txt -o go_by_genus.biom
```

The parameter `-f 2` instructs the program to collapse by the second field of each feature (e.g., `K00133` of `Escherichia|K00133`). See [here](collapse.md#stratification) for details.
