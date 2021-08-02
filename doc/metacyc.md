# Working with MetaCyc

**MetaCyc** (https://metacyc.org/) ([Caspi et al., 2020](https://academic.oup.com/nar/article/48/D1/D445/5581728)) is a metabolic pathway database that has been widely used in genomic, metagenomic and metabolomic studies. It provides a hierarchical classification system, including genes, proteins, reactions, compounds, pathways and more.


## Contents

- [Mapping files](#mapping-files)
- [Protein profiling](#protein-profiling)
- [Genes, reactions, compounds and pathways](#genes,-reactions,-compounds-and-pathways)
- [Pathway coverage](#pathway-coverage)


## Mapping files

We mapped all ORFs from the [WoL](wol.md) reference genome database to the reference protein sequences in MetaCyc release 23.0. We provide this mapping file, as well as Woltka-compatible mapping and annotation files representing the higher levels in the MetaCyc classification system corresponding to the mapped WoL ORFs. These files are publicly available under the `function/metacyc/` directory of the WoL data release ([see details](wol.md)).

We also included a UniRef-to-MetaCyc mapping file, extracted from the [UniProt](https://www.uniprot.org/downloads) data release and subsetted to WoL. It contains less entries though.

The original, full MetaCyc database is available for download from the [official website](https://metacyc.org/). It includes the reference protein sequences with which one can perform custom alignments. We provide a Python script: [**metacyc_build.py**](https://github.com/qiyunzhu/utils/blob/main/metacyc_build.py) to reformat a local copy of the MetaCyc database into simple mapping files which can be parsed by Woltka.


## Protein profiling

The following command utilizes Woltka's [ordinal classification](ordinal.md) function to classify sequences to **proteins** -- the entry point of the MetaCyc hierarchies.

```bash
woltka classify \
  --input  input_dir \
  --coords coords.txt.xz \
  --map    metacyc/protein.map.xz \
  --names  metacyc/protein_name.txt \
  --rank   protein \
  --output protein.biom
```

* Note: The `--names` parameter is optional.
* Note: Replacing `.biom` with `.tsv` can generate TSV output.

Alternatively, one can split this command into two:

```bash
woltka classify -i input_dir -c coords.txt.xz -o orf.biom
woltka tools collapse -i orf.biom -m metacyc/protein.map.xz -n metacyc/protein_name.txt -o protein.biom
```


## Genes, reactions, compounds and pathways

The MetaCyc hierarcies are as follows:

```
               v
       go < protein > gene > pathway
               v
regulation < enzrxn
               v
       ec < reaction > compound (left / right) > type
               v
     type < pathway > taxonomic range
               v
         super pathway
               v
              type
```

All transitions are enabled using Woltka's [collapse](collapse.md) command with individual mapping files. For example, one can generate profiles along the following cascade:

```
protein - enzrxn - reaction - pathway - super pathway - type
```

Using the following commands:

```bash
# protein to enzrxn (enzymatic reaction):
woltka tools collapse -i protein.biom -m metacyc/protein-to-enzrxn.txt -n metacyc/enzrxn_name.txt -o enzrxn.biom
# enzrxn to reaction:
woltka tools collapse -i enzrxn.biom -m metacyc/enzrxn-to-reaction.txt -n metacyc/reaction_name.txt -o reaction.biom
# reaction to pathway:
woltka tools collapse -i reaction.biom -m metacyc/reaction-to-pathway.txt -n metacyc/pathway_name.txt -o pathway.biom
# pathway to super pathway:
woltka tools collapse -i pathway.biom -m metacyc/pathway-to-super_pathway.txt -n metacyc/pathway_name.txt -o super_pathway.biom
# super pathway (or pathway) to pathway type:
woltka tools collapse -i super_pathway.biom -m metacyc/pathway_type.txt -n metacyc/all_class_name.txt -o pathway_type.biom
```

The collapse command supports **many-to-many** mapping. For example, if one reaction is found in three pathways, each pathway will be counted **once**. In some instances (e.g., to retain compositionality of the profile), one may consider adding the `--frac` flag, which will instruct the program to count each pathway 1 / 3 times ([see details](collapse.md)).


## Pathway coverage

It is usually important to assess how **completed**, in addition to how abundant a metabolic pathway is in a sample. This is because a pathway can only function when all components are present. Woltka's [coverage](coverage.md) command provides a solution ([see details](coverage.md)).

Among the mapping files there are two files recording the composition of metabolic pathways in terms of **genes** and **reactions**. We recommend using the reaction mapping because there are duplicated gene IDs.

```bash
woltka tools coverage -i reaction.biom -m pathway-to-reaction_list.txt -o pathway_coverage.biom
```

The output file is a **coverage** table in which every cell value represents the percentage of member reactions of a particular pathway present in a particular sample.


## Stratifying by taxonomy

Woltka allows overlaying two or more classification schemes on the same profile using the [stratification](stratify.md) function. With this function, one can identify functional genes and metabolic pathways that belong to individual taxonomic groups.

First, perform taxonomic classification at the desired rank (e.g., genus) and generate read-to-genus maps (using parameter `--outmap` or `-u`).

```bash
woltka classify -i input_dir --lineage lineages.txt -r genus -o genus.biom -u map_dir
```

Second, perform functional classification. This command is identical to the first command in this document, except for the addition of `--stratify` or `-t` parameter pointing to the genus maps, which will be incorporated into the functional classes ([see details](stratify.md)).

```bash
woltka classify -i input_dir -c coords.txt.xz -m metacyc/protein.map.xz -r protein -t map_dir -o protein.biom
```
