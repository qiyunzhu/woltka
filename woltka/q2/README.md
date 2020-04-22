# q2-woltka

**q2-woltka** is a QIIME 2 plugin for Woltka, the metagenome data analysis tool.

We recommend that you use Woltka with the WoL reference phylogeny ([**tree.qza**](https://biocore.github.io/wol/data/trees/tree.qza)) and taxonomy ([**taxonomy.qza**](https://biocore.github.io/wol/data/taxonomy/ncbi/taxonomy.qza)). However we note that Woltka is highly flexible and it supports custom classification systems.

## Installation

Requires: QIIME 2 2019.1 or above.

First, install the standalone Woltka program:

```bash
pip install git+https://github.com/qiyunzhu/woltka.git
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

Woltka ships with small test datasets under this directory:

```
<program_dir>/woltka/tests/data
```

The following command can generate a gOTU table based on per-sample sequence alignment files under the same directory:

```bash
qiime woltka gotu --p-input-path align/bowtie2 --o-table table.qza
```

In a real scenario, the workflow for shotgun metagenome sequence alignment is available at the Qiita server. One may process "per_sample_FASTQ" files with command: "**Shogun v1.0.7**" and parameter set: "**wol_bowtie2**". This will launch the alignment job. Once done, one may download the alignment file `alignment.bowtie2.sam.xz`, then generate a gOTU table with:

```bash
qiime woltka gotu --p-input-path alignment.bowtie2.sam.xz --o-table table.qza
```

The output file `table.qza` contains a BIOM table (QIIME artifact type: `FeatureTable[Frequency]`) in which columns are samples and rows are gOTUs. It can then be analyzed using the classical QIIME 2 workflow.
