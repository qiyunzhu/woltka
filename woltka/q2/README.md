# q2-woltka

**q2-woltka** is a QIIME 2 plugin for Woltka. It retains the core functionalities of Woltka, while its command-line interface was modified to be consistent with the design goal and other plugins of QIIME 2.

We recommend that you use Woltka with the WoL reference phylogeny ([tree.qza](https://biocore.github.io/wol/data/trees/tree.qza)) and taxonomy ([taxonomy.qza](https://biocore.github.io/wol/data/taxonomy/ncbi/taxonomy.qza)). However we note that Woltka is highly flexible and it supports custom classification systems.


## Contents

- [Installation](#installation)
- [Example usage](#example-usage)
- [Data importing](#data-importing)
- [Resulting profile](#resulting-profile)
- [Command-line interface](#command-line-interface)
- [Difference from Woltka](#difference-from-woltka)


## Installation

Requirement: **QIIME 2 2019.1** or above.

First, install the standalone Woltka program:

```bash
pip install woltka
```

Second, update QIIME 2 cache and the plugin will be automatically loaded into QIIME 2.

```bash
qiime dev refresh-cache
```

Then one can launch the program by executing:

```bash
qiime woltka
```

## Example usage

Multiple small test files can be found under [tests/data](tests/data). To access them, [download](https://github.com/qiyunzhu/woltka/archive/master.zip) this GitHub repo, unzip, and navigate to this directory.

[**Note**] These files are QIIME 2 artifacts (*.qza). They are containers of actual data files. See [below](#data-importing) for how to generate these files.


Generate an OGU (operational genomic unit) table:

```bash
qiime woltka classify \
  --i-alignment align_sam.qza \
  --o-classified-table table.qza
```

Classify a SAM alignment file at the genus rank using a QIIME 2 taxonomy (lineage strings):

```bash
qiime woltka classify \
  --i-alignment align_sam.qza \
  --p-target-rank genus \
  --i-reference-taxonomy taxonomy.qza \
  --o-classified-table table.qza
```

Perform free-rank classification on a BLAST 6 output file using an NCBI nodes file:

```bash
qiime woltka classify \
  --i-alignment align_b6o.qza \
  --p-target-rank free \
  --i-taxon-map g_to_taxid.qza \
  --i-reference-nodes nodes.qza \
  --o-classified-table table.qza
```

Perform phylogenetic classification on a simple mapping file using a QIIME 2 phylogeny (Newick tree):

```bash
qiime woltka classify \
  --i-alignment align_map.qza \
  --p-target-rank free \
  --p-subject-is-okay \
  --i-reference-tree tree.qza \
  --o-classified-table table.qza
```

Perform coordinates-based functional classification at the pathway rank using an nodes-formatted MetaCyc hierarchy:

```bash
qiime woltka classify \
  --i-alignment align_sam.qza \
  --p-target-rank pathway \
  --i-gene-coordinates coords.qza \
  --i-taxon-map gene_to_uniref.qza \
  --i-reference-nodes metacyc_nodes.qza \
  --o-classified-table table.qza
```


## Data importing

QIIME 2 [artifacts](https://docs.qiime2.org/2020.6/concepts/?highlight=artifact#data-files-qiime-2-artifacts) are containers of data files. To use **q2-woltka**, one first needs to [import](https://docs.qiime2.org/2020.6/tutorials/importing/) data files into artifacts. Example:

```bash
qiime tools import --type FeatureData[SeqAlnMap] --input-path alignment.sam --output-path alignment.qza
```

In this example, the [semantic type](https://docs.qiime2.org/2020.6/concepts/?highlight=artifact#semantic-types) of the data file is `FeatureData[SeqAlnMap]`, which means the SAM format.

q2-woltka defines five new semantic types:

- `FeatureData[SeqAlnMap]`: [Sequence alignment map (SAM) format](../../doc/input.md#file-formats).
- `FeatureData[BLAST6Out]`: [BLAST tabular output format 6](../../doc/input.md#file-formats).
- `FeatureData[SimpleMap]`: [Simple mapping format](../../doc/input.md#file-formats) (and [here](../../doc/hierarchy.md#supported-hierarchy-files)).
- `FeatureData[NCBINodes]`: [NCBI nodes format](../../doc/hierarchy.md#supported-hierarchy-files).
- `GeneCoordinates`: [Gene coordinate mapping format](../../doc/ordinal.md#gene-coordinates).


## Downstream analyses

The resulting file (e.g., `table.qza`) is a **feature table**. Its semantic type is `FeatureTable[Frequency]`. This table can be directly analyzed using downstream QIIME plugins. For example:

```bash
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny tree.qza \
  --i-table table.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file metadata.tsv \
  --output-dir .
```


## Command-line interface

The following table lists all parameters under `q2-woltka classify` and their corresponding [Woltka commands](../../doc/cli.md). Note that not all Woltka commands are implemented in q2-woltka, and the behaviors of some commands are altered (detailed [below](#difference-from-woltka)).

### Inputs

q2-woltka option | Woltka option | Description
--- | --- | ---
`--i-alignment` | `-i\|--input` | Multiplexed sequence alignment map to be classified. Can accept SAM, BLAST6 or simple map format (required).
`--i-reference-taxonomy` | `--lineage` | Reference taxonomic lineage strings.
`--i-reference-tree` | `--newick` | Reference phylogenetic tree.
`--i-reference-nodes` | `--nodes` | Reference taxonomic nodes.
`--i-taxon-map` | `-m\|--map` | Mapping of subject IDs to taxon IDs.
`--i-gene-coordinates` | `-c\|--coords` | Reference taxonomic lineage strings.

### Parameters

q2-woltka option | Woltka option | Description
--- | --- | ---
`--p-target-rank` | `-r\|--rank` | Classify sequences at this rank. Enter "none" to directly report subjects; enter "free" for free-rank classification (required).
`--p-trim-subject` | `--trim-sub` | Trim subject IDs at the last underscore.
`--p-overlap-threshold` | `--overlap` | Read/gene overlapping percentage threshold.
`--p-unique-assignment` | `--uniq` | One sequence can only be assigned to one classification unit, or remain unassigned if there is ambiguity. Otherwise, all candidate units are reported and their counts are normalized.
`--p-majority-threshold` | `--major` | In given-rank classification, use majority rule at this percentage threshold to determine assignment when there are multiple candidates.
`--p-above-given-rank` | `--above` | In given-rank classification, allow assigning a sequence to a higher rank if it cannot be assigned to the current rank.
`--p-subject-is-okay` | `--subok` | In free-rank classification, allow assigning a sequence to its direct subject, if applicable, before going up in hierarchy.
`--p-report-unassigned` | `--unassigned` | Report Frequency of unassigned sequences (will be marked as "Unassigned").

### Outputs

q2-woltka option | Woltka option | Description
--- | --- | ---
`--o-classified-table` | `-o\|--output` | The resulting table of frequencies of classification units (required).


## Difference from Woltka

**q2-woltka** includes core features of the standalone Woltka. The results are identical, if using corresponding settings. But some controls are not provided, while the interface and the behaviors of some parameters are not identical. The differences are listed below.

- The input alignment must be **multiplexed**, with sample ID and an underscore (`_`) prefixing every read ID.

- Only one profile is generated per run. The profile must be in BIOM format. The parameter `--p-target-rank` (`-r`) only accepts one value and cannot be omitted.

- Only one classification system is allowed per run. The user can choose from:
  1. "**taxonomy**" (Greengenes-style lineage strings, i.e., QIIME 2's built-in semantic type `FeatureData[Taxonomy]`),
  2. "**tree**" (Newick format tree, i.e., QIIME 2's built-in semantic type `Phylogeny[Rooted]`]), and
  3. "**nodes**" (NCBI-style nodes, i.e., the new semantic type `FeatureData[NCBINodes]`).

- Optionally, a subject-to-taxon mapping can be added. But only one map is allowed and this map does not bear rank information.

[**Note**] In order to pursue the flexibility of the standalone Woltka, one needs to customize a **nodes** file to represent the classification system.

- Taxon names and ranks are not appended to the output BIOM table (not allowed in QIIME 2).

- Currently does not have the feature map writing function.

- Currently does not have the stratification function.
