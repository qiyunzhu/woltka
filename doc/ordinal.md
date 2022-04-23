# Coordinates-based classification

Woltka combines the two fundamental analyses in metagenomics: Structural profiling (mapping reads to **genomes**) and functional profiling (mapping reads to functional **genes**) into one run. This saves compute, ensures consistency, and allows for stratification which is essential for understanding functional diversity across the microbial community.

This is achieved by an efficient algorithm implemented in Woltka, which matches read alignments and annotated genes based on their **coordinates** on the host genome. The algorithm matches reads and genes in an [ordinal](https://en.wikipedia.org/wiki/Ordinal_number) scale, with the efficiency of _O_(_n_+_m_) (_n_ and _m_ are numbers of reads and genes), therefore suitable for handling large microbiome datasets and databases.

## Contents

- [Input alignments](#input-alignments)
- [Gene coordinates](#gene-coordinates)
- [Classification system](#classification-system)
- [Matching threshold](#natching-threshold)


## Input alignments

The alignment files must contains the information of alignment coordinates. Both SAM and BLAST formats (detailed [here](input.md)) have this information. They can just be supplied in the regular way.


## Gene coordinates

One needs to provide Woltka with a table of gene coordinates on host genomes. Supported formats include:

Native NCBI annotations and accessions:

```
## GCF_000005825.2
# NC_013791.2
WP_012957018.1  816 2168
WP_012957019.1  2348    3490
WP_012957020.1  3744    3959
WP_012957021.1  3971    5086
...
```

WoL reannotations (available for [download](https://biocore.github.io/wol/)), in which nucleotides of each genome are concatenated, and a numeric index is used to identify each gene. Woltka will automatically prefix these numbers with the genome ID, such as `G000006745_1`.

```
>G000006745
1       806     372
2       2177    816
3       3896    2271
4       4446    4123
5       4629    4492
```

Compressed files are supported.


## Classification system

The same as the regular ones detailed [here](hierarchy.md), but subject IDs are gene IDs instead of genomes. For examples, `WP_012957018.1` or `G000006745_1`.


## Matching threshold

Woltka considers a read-to-gene match if their **overlapping** region reaches a percentage threshold of the entire alignment. This is controled by the parameter `--overlap`. The default value is 80 (%).


## Example command

With the coordinates file, one can streamline the read-to-gene matching step into a Woltka protocol. Here is an example for functional profiling:

```bash
woltka classify \
  --coords coords.txt \
  --overlap 95 \
  --map gene2function.txt \
  --map function2pathway.txt \
  --rank function,pathway \
  ...
```
