## Coordinates-based functional analyses

Woltka combines the two fundamental analyses in metagenomics: taxonomic profiling (mapping reads to genomes) and functional profiling (mapping reads to functional genes) into one run. This saves compute, ensures consistency, and allows for stratification which is essential for understanding functional diversity across the microbial community.

This is achieved by an efficient algorithm implemented in Woltka, which matches read alignments and annotated genes based on their coordinates on the host genome.

The coordinates of read-to-genome alignments are provided in the alignment files. One needs to provide Woltka with a table of gene coordinates. The format is like:

WoL reannotations (available for [download](https://biocore.github.io/wol/)):

```
>G000006745
1       806     372
2       2177    816
3       3896    2271
4       4446    4123
5       4629    4492
```

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

Again, compressed files are supported.

With the coordinates file, one can streamline the read-to-gene matching step into a Woltka protocol. Here is an example for functional profiling:

```bash
woltka classify \
  --coords coords.txt \
  --map gene2function.txt \
  --map function2pathway.txt \
  ...
```
