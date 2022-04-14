# Profile normalization

By default, the cell values in a feature table (profile) are **counts** (**frequencies**) of each feature in each sample. This behavior maximizes the flexibility of downstream analyses. Meanwhile, Woltka provides multiple functions for normalizing cell values to alternative units and scales:

- `--sizes`: Divide counts by subject sizes (e.g., by genome or gene lengths).
- `--frac`: Convert counts into fractions (relative abundances).
- `--scale`: Scale up all values by a constant factor (e.g., "1k", "1M").
- `--digits`: Keep certain number of digits after the decimal point.


## Contents

- [During or after classification](#during-or-after-classification)
- [Normalization by subject size](#normalization-by-subject-size)
  - [Sequence and taxonomic abundances](#sequence-and-taxonomic-abundances)
  - [Why real-time normalization matters](#Why-real-time-normalization-matters)
  - [Abundance of functional genes](#abundance-of-functional-genes)
- [Relative abundance](#relative-abundance)


## During or after classification

Woltka has two modes of normalization, each providing all four functions mentioned above.

- **First**, The normalization functions can be called as part of the main classification workflow (`woltka classify`), and directly output normalized profiles.

- **Alternatively**, one can perform normalization on profiles that are already generated (`woltka tools normalize`). This saves the need for rerunning the lengthy classification process.

There is a major advantage of normalization by subject size **during** classification, as detailed below. Other than that, the two modes produce mutually identical results. See below for some examples.


## Normalization by subject size

The number of reads assigned to each classification unit can be normalized against a **subject-specific** property, which is usually size (length), but can also be other metrics depending on the specific research.

During classification (as part of main workflow):

```bash
woltka classify --sizes size.map ...
```

Or post classification (on existing profiles):

```bash
woltka tools normalize --sizes size.map ...
```

### Sequence and taxonomic abundances

A common usage of this function is to convert **sequence abundance** to **taxonomic abundance**.

Woltka's taxonomic classification protocol reports sequence abundance (i.e., the _number of sequencing reads_ assigned to each taxonomic group). However, in some scenarios, the taxonomic abundance (the _number of genomes_ per taxonomic group found in the sample) is more relevant to the underlying biology.

- **Note**: It is important to understand the difference between these two abundances. In [Sun et al. (2021)](https://www.nature.com/articles/s41592-021-01141-3), our collaborators and us demonstrated that confusing the two notions will lead to highly biased profiling results and biological conclusions.

To perform this convertion, one needs a genome ID to genome length (bp) mapping file. The [WoL data release](wol.md) provides this file as `length.map`.

It is more convenient to scale up resulting values by a constant, such as one million (i.e., **copies per million, or CPM**). This can be achieved using the `--scale` parameter.

Examples (post classification):

```bash
woltka tools normalize -i ogu.biom -m length.map --scale 1M -o ogu.cpm.biom
```

During classification:

```bash
woltka classify \
  --input  input.sam \
  --map    taxonomy/taxid.map \
  --nodes  taxonomy/nodes.dmp \
  --names  taxonomy/names.dmp \
  --rank   none,free,species,genus,phylum \
  --sizes  length.map \
  --scale  1M \
  --output outdir
```

### Why real-time normalization matters

What's the difference between the two modes? As we can see, if normalization is performed during classification, one can generate multiple rank-specific profiles all at once.

But the benefit is not just convenience. In order to normalize an existing profile, one has to provide a feature-to-size mapping. It probably isn't hard to find the lengths of reference genomes (bp), however, what about the "size" of taxonomic units?

One may think of workarounds such as averaging the lengths of genomes under each taxonomic unit. This may work, but has significant limitations, as genome lengths vary greatly, especially with highly elastic microbial lineages (such as _Escherichia coli_) and higher taxonomic ranks, taxon sampling bias in the database, and assembly artifacts in the reference genomes.

In contrast, Woltka directly normalizes counts of reference genomes **before** they are classified into highly taxonomic units. The result reflects the **exact** taxonomic abundance: the quantity of genomes under each taxonomic group.

For example, in a metagenome, 100 reads are assigned to genome A (2.0 Mb); 200 reads are assigned to genome B (5.0 Mb). Both genomes A and B are classified as species X. Then the normalized abundance of species X is 100 / 2.0 + 200 / 5.0 = 90 CPM.

### Abundance of functional genes

This function is also useful in normalizing a functional profile by the lengths of genes, which is a common practice in the field of transcriptomics, and has seen popularity in metagenomics (e.g., [HUMAnN](https://github.com/biobakery/humann)).

Woltka provides a convenient shortcut: The gene coordinates file (see [details](ordinal.md)) can be directly used in replacement of an explicit gene length mapping file.

```bash
woltka classify \
  --input  input_dir \
  --coords coords.txt.xz \
  --map    metacyc/protein.map.xz \
  --names  metacyc/protein_name.txt \
  --rank   protein \
  --sizes   . \
  --scale  1k \
  --digits 3 \
  --output protein.tsv
```

The output values are in the unit of **reads per kilobase, or RPK**. They reflect the quantity of functional genes found in the sample.


### Relative abundance

One can convert counts into **fractions** of total counts of each sample, i.e., **relative abundances**, which are useful in many applications.

During classification:

```bash
woltka classify ... --frac
```

The output values are fractions (like "0.05"). One can add `--scale 100` to convert them to **percentages** (like "5" (%)).

On an existing profile:

```bash
woltka tools normalize ... (without --sizes, no need for --frac)
```

**Note**: These values are the fractions of reads assigned to each classification units versus all reads that are **assigned**. If you add `--unassigned`, the values become the fractions of reads versus all reads that are **aligned**. Woltka cannot calculate the fractions of reads versus the **original** sequencing data, since it processes alignment files instead of raw FastQ files. However, it isn't hard to do this calculation manually if you know the sequencing depth information.
