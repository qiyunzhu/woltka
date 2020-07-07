# Types of classifications

Woltka grants users the flexibility to control the classification criteria, in addition to the flexibility in the classification system [itself](hierarchy.md). It is usually preferrable to run Woltka using alternative parameters, and see which ones work the best



## Contents

- [Target rank (or no rank)](#target-rank-or-no-rank)
- [How Woltka handles hierarchies](#how-woltka-handles-hierarchies)
- [Feature name dictionary](#feature-name-dictionary)


## Target rank (or no rank)

Woltka features the following modes of classification, as controlled by the `--rank` parameter:

1\. No classification (`--rank none`):

2\. Free-rank classification (`--rank free`)

3\. Given-rank classification (`--rank <name>`). The name can be `species`, `genus`, `family`..., or `K`, `M`, `R`, `map`..., or `ATC4`, `ATC3`..., or whatever, as long as the classification [hierarchies](hierarchy.md#supported-hierarchy-files) you supplied have them.

Multiple ranks can be specified simultaneously, delimited by comma (e.g., `none,free,phylum,genus,species`), in which case Woltka will generate one [profile](output.md) for each rank.

Woltka features a highly flexible hierarchical classification system. It is represented by a **tree** structure instead of a fixed number of levels (e.g., the eight standard taxonomic ranks). In another word, it is **rank-free**.


## Lowest common ancestor


## Majority rule


## Unassigned