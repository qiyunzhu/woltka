# Classification methods

Woltka grants users the flexibility to control the classification criteria, in addition to the flexibility in the classification system [itself](hierarchy.md). This maximizes the potential of exploring microbiome big data in various angles, and enables novel applications.


## Contents

- [Target rank (or no rank)](#target-rank-or-no-rank)
- [Ambiguous assignment](#ambiguous-assignment)
- ["Unassigned" sequences](#"unassigned"-sequences)


## Target rank (or no rank)

Woltka features the following modes of classification, as controlled by the `--rank` (or `-r`) parameter:

### 1\. **No** classification (`--rank none`)

Simply report subject IDs. A classification system is not required in this analysis. This mode has the highest granularity and is useful in e.g. the [gOTU](gotu.md) analysis.

### 2\. **Free-rank** classification (`--rank free`)

Find the best classification unit in the entire hierarchy in describing the query sequence, without forcing it to a particular rank.

This mode uses a lowest common ancestor (LCA) algorithm designed for a tree structure without fixed ranks to tackle ambiguous assignments. It chooses the lowest unit when possible, and go higher in the hierarchy if necessary.

### 3\. **Given-rank** classification (`--rank <name>`)

Choose a classification unit at the given rank to describe the query sequence. This is the closest to the conventional notion of "taxonomic classification".

The rank can be `species`, `genus`, `family`..., or `K`, `M`, `R`, `map`..., or `ATC4`, `ATC3`..., or whatever, as long as the classification [hierarchies](hierarchy.md#supported-hierarchy-files) you supplied have them.

### 4\. **Given-rank-and-above** classification (add flag `--above` in addition to `--rank <name>`)

Attempt to classify at the given rank, and when this is not possible, go up in hierarchy until a proper unit is reached. The same LCA algorithm is used in this process.

The difference from free-rank classification (2) is that, it discards all units below the given rank, even if some of them may describe the query sequence.

### Multiple ranks

Multiple ranks can be specified simultaneously, delimited by comma (e.g., `--rank none,free,phylum,genus,species`), in which case Woltka will generate one [profile](output.md) for each rank. This is significantly faster than running Woltka multiple times on individual ranks.


## Ambiguous assignment

In many cases a query sequence has matches in multiple reference sequences, and those sequences may have have different assignments at the rank you instruct Woltka to classify at. Woltka deals with this situation using your choice of the following mechanisms.

### 1\. (Default mode) keep them all, and normalize

In the resulting profile, each subject receives 1 / _k_ count, where _k_ is the total number of subjects of the current query.

For example, sequence A was aligned to five genomes: two under genus [_Escherichia_](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=561&lvl=3&lin=f&keep=1&srchmode=1&unlock), one under each of genera [_Salmonella_](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=590&lvl=3&lin=f&keep=1&srchmode=1&unlock), [_Klebsiella_](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=570&lvl=3&lin=f&keep=1&srchmode=1&unlock), and [_Enterobacter_](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=547&lvl=3&lin=f&keep=1&srchmode=1&unlock). So at **genus** level, _Escherichia_ receives 2/5 count, and each of the other three receives 1/5 each.

In the read-to-feature maps (`--outmap`), a multi-assignment will be reported as:

```
A <tab> Escherichia:2 <tab> Salmonella:1 <tab> Klebsiella:1 <tab> Enterobacter:1
```

### 2\. Unique assignment (`--uniq`)

It applies to none (1) or given-rank (3) classifications. With this flag, ambiguous assignments will be considered as **unassigned**.

In the above case, this query sequence won't receive any genus-level assignment due to the ambiguity.

### 3\. Lowest common ancestor ([LCA](https://en.wikipedia.org/wiki/Lowest_common_ancestor))

This is how free-rank (2) and given-rank-and-above (4) classifications work. No parameter is needed.

In the above case, the query sequence will be assigned to family [Enterobacteriaceae](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=543&lvl=3&lin=f&keep=1&srchmode=1&unlock), because all four genera belong to this family.

### 4\. Majority rule (`--major <%>`)

At a given rank (3), as long as the dominant unit reaches the given percentage threshold of all subjects, it will be considered as the right target.

In the above case, at the family level, all four genera (five genomes) belong to Enterobacteriaceae (100%), so there is no doubt that it can be assigned to Enterobacteriaceae. Now we replace one genome with one under genus [_Pseudomonas_](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=286&lvl=3&lin=f&keep=1&srchmode=1&unlock) (family [Pseudomonadaceae](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=135621&lvl=3&lin=f&keep=1&srchmode=1&unlock)). Therefore, the proportion of Enterobacteriaceae becomes 4/5 = 80%.

The query sequence will be assigned to Enterobacteriaceae if you specify `--major 80` or above, or unassigned if this parameter is omitted or below 80.

Note: Majority rule overrides the `--above` flag. Currently, Woltka cannot combine LCA and majority rule due to the rank-free nature of the classification system.


## "Unassigned" sequences

With flag `--unassigned`, Woltka reports unassigned sequences in the profile and the feature map. They will be marked as "Unassigned".

A sequence is deemed unassigned because of one of the following reasons:

1. The subject(s) is **not found** in the classification system. For flexibility, Woltka does NOT consider this as a conflict between data and database, and it does NOT halt the program and warn the user. Instead, it is treated as unassigned.

2. The LCA (see above) of subjects is the **root**. A [tree-structured](hierarchy.md) classification hierarchy always have a root, and no matter how diverse subjects are, they always coalesce to the root eventually. But reporting "root" as an assignment is meaningless. So it will be considered as unassigned.

3. In unique-assignment mode, assignments of subjects are not unique (see above).

4. In majority-rule mode, none of the candidate units reaches the threshold (see above).

[**Important**] The "unassigned" part represent query sequences that were **aligned** to one or more **subjects**, but Woltka cannot find a suitable assignment based on those subjects. Therefore, assigned + unassigned is NOT the entire sample, but only the part of sample that are found in the alignment.

In [ordinal mapping](ordinal.md), "subjects" are **genes** instead of genomes, therefore despite that some query sequences can be aligned to one or more genomes, they can still be excluded from the "unassigned" part if their coordinates do not match any gene.
