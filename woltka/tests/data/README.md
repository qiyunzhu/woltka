# Test Datasets and Commands

The test datasets, including reference genome sequences, taxonomy, and input sequencing data are based on a a collection of 107 NCBI-defined ["reference"](https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes) bacterial genomes.

## Align

**Sequence alignment** results.

Here "alignment" refers to the operation of aligning short DNA sequences (**reads**) against reference **genome** sequences.

The query sequences are five samples (S01 to S05) of 150 bp paired-end reads simulated using [CAMISIM](https://github.com/CAMI-challenge/CAMISIM/). The ground-truth mapping of reads against original genomes and locations are provided in `truth`.

Five aligners were used: [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Centrifuge](https://ccb.jhu.edu/software/centrifuge/), [DIAMOND](http://diamondsearch.org/index.php), and [BURST](https://github.com/knights-lab/BURST). In particular, `bowtie2` used the Bowtie2 program with its default setting, whereas `bt2sho` used the Bowtie2 parameters recommended by SHOGUN. The mappings between reads and genomes are provided by each alignment file. `diamond` mapped reads to reference genes instead of genomes.

In addition, `burst/split` is the mapping against **genes** annotated from the reference genomes (not genomes themselves).

## Taxonomy

**Reference taxonomy** of the 107 bacterial genomes. The taxonomy is provided in multiple formats:

- `taxid.map`: Genome ID to NCBI taxonomy ID mapping file.
- `taxdump`: NCBI taxdump-style database files (`nodes.dmp` and `names.dmp`).
- `lineages.txt`: Greengenes-style lineage strings.
- `rank_names.tsv` and `rank_tids.tsv`: Taxon name and TaxID at each of the seven standard taxonomic ranks for each genome.

In addition, `nucl/` contains the mappings from nucleotide sequence accessions (instead of their host genomes) to taxonomy.

## Phylogeny

- `tree.nwk`: Maximum likelihood phylogenetic tree of the 107 bacterial genomes built using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) and midpoint-rooted.

## Function

**Functional annotations** of genes ([Prodigal](https://github.com/hyattpd/Prodigal)-predicted ORFs) on the 107 genomes.

- `coords.txt`: Coordinates of genes on their host genomes
- `uniref.map`: Mapping of genes to [UniRef](https://www.uniprot.org/help/uniref) entries.
- `go/`: Mapping of UniRef entries to [GO](http://geneontology.org/docs/ontology-documentation/) terms.
- `nucl/`: See above.

## Output

Woltka output files generated from the provided alignments and classification systems. The commands for generating those files are:

`bowtie2.ogu.tsv`: Generate an [OGU table](../../../doc/ogu.md) (just assign, don't classify).

```bash
woltka classify -i align/bowtie2 -o bowtie2.ogu.tsv
```

`bowtie2.free.tsv`: [Free-rank](../../../doc/classify.md#2-free-rank-classification---rank-free) classification

```bash
woltka classify \
  --input align/bowtie2 \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --rank free \
  --output bowtie2.free.tsv
```

`bowtie2.free.1p.tsv`: [Filter](../../../doc/filter.md) out features with <1% per-sample abundance.

```bash
woltka tools filter \
  --input output/bowtie2.free.tsv \
  --min-percent 1 \
  --output bowtie2.free.1p.tsv
```

`blastn.species.tsv`: Taxonomic profiling of a [multiplexed](../../../doc/input.md#demultiplexing) alignment file using a Greengenes-style [lineage](../../../doc/hierarchy.md#3-lineage-strings---lineage) file.

```bash
woltka classify \
  --input align/blastn/mux.b6o.xz \
  --lineage taxonomy/lineages.txt \
  --rank species \
  --output blastn.species.tsv
```

`blastn.species.biom`: Same as above, replacing `.tsv` with `.biom`.

`burst.genus.tsv`:

```bash
woltka classify \
  --input align/burst \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank genus \
  --name-as-id \
  --output burst.genus.tsv
```

`burst.genus.map/`: Write [read-to-genus mapping](../../../doc/output.md#output-read-maps) files.

Same as above, adding:

```bash
--outmap burst.genus.map
```

`bt2sho.phylo.tsv`: Classify sequences using a [phylogenetic tree](../../../doc/hierarchy.md#2-newick-tree---newick). Resulting features are tips and internal nodes of the tree.

```bash
woltka classify \
  --input align/bt2sho \
  --newick tree.nwk \
  --rank free \
  --name-as-id \
  --output bt2sho.phylo.tsv
```

`bt2sho.order.cpm.tsv`: Feature frequencies are [normalized by genome length](../../../doc/normalize.md#normalization-by-subject-size), transformed to the unit of CPM (counts per million), rounded to 3 digits after the decimal point.

```bash
woltka classify \
  --input align/bt2sho \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank order \
  --sizes taxonomy/length.map \
  --scale 1M \
  --digits 3 \
  --output bt2sho.order.cpm.tsv
```

`burst.process.tsv`: Functional profiling to GO biological processes.

```bash
woltka classify \
  --input align/burst \
  --coords function/coords.txt.xz \
  --map function/uniref/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --rank process \
  --output burst.process.tsv
```

`burst.genus.process.tsv`: [Stratify](../../../doc/stratify.md) functional units by taxonomic units.

Same as above, adding:

```bash
--stratify burst.genus.map
```

`split.genus.tsv`: Input are read-to-gene alignments. Translate gene IDs to host chromosome IDs (by trimming at `_`), then perform taxonomic profiling.

```bash
woltka classify \
  --input align/burst/split \
  --trim-sub _ \
  --map taxonomy/nucl/nucl2tid.txt \
  --nodes taxonomy/nodes.dmp \
  --names taxonomy/names.dmp \
  --rank genus \
  --name-as-id \
  --output split.genus.tsv
```

`split.process.tsv`: Input are read-to-gene alignments. Perform functional profiling.

```bash
woltka classify \
  --input align/burst/split \
  --map function/nucl/uniref.map.xz \
  --map function/go/process.tsv.xz \
  --rank process \
  --output split.process.tsv
```

`merged.process.tsv`: [Merge](../../../doc/merge.md) two profiles.

```bash
woltka tools merge \
  --input output/burst.process.tsv \
  --input output/split.process.tsv \
  --output merged.process.tsv
```

`diamond.free.tsv`: Input are read-to-protein alignments. Perform taxonomic profiling.

```bash
woltka classify \
  --input align/diamond \
  --trim-sub _ \
  --map taxonomy/taxid.map \
  --nodes taxonomy/nodes.dmp \
  --rank free \
  --output diamond.free.tsv
```

`diamond.function.tsv`: Input are read-to-protein alignments. Perform functional profiling.

```bash
woltka classify \
  --input align/diamond \
  --map function/uniref/uniref.map.xz \
  --map function/go/function.tsv.xz \
  --rank function \
  --output diamond.function.tsv
```

`bt2sho.component.rpk.tsv`: Perform functional profiling, [normalize](../../../doc/normalize.md) gene frequencies by gene lengths, then transform to the unit of RPK (reads per kilobases).

```bash
woltka classify \
  --input align/bt2sho \
  --coords function/coords.txt.xz \
  --map function/uniref/uniref.map.xz \
  --map function/go/component.tsv.xz \
  --rank component \
  --sizes . \
  --scale 1k \
  --digits 3 \
  --output bt2sho.component.rpk.tsv
```

`truth.gene.tsv`: Input are ground truth. Perform simple [read-gene matching](../../../doc/ordinal.md).

```bash
woltka classify \
  --input align/truth/b6o \
  --coords function/nucl/coords.txt.xz \
  --output truth.gene.tsv
```

`truth.uniref.tsv`: Perform functional profiling by UniRef entries.

```bash
woltka classify \
  --input align/truth/b6o \
  --coords function/nucl/coords.txt.xz \
  --map function/nucl/uniref.map.xz \
  --names function/uniref/uniref.names.xz \
  --rank uniref \
  --output truth.uniref.tsv
```

Or, instead of doing the one-step classification workflow, [collapse](../../../doc/collapse.md) an existing per-gene profile into a UniRef profile.

```bash
woltka tools collapse \
  --input output/truth.gene.tsv \
  --map function/nucl/uniref.map.xz \
  --names function/uniref/uniref.names.xz \
  --output truth.uniref.tsv
```

`truth.goslim.tsv`: Further collapse the UniRef profile into a [GO Slim](https://www.ebi.ac.uk/QuickGO/help/slims) profile. The program can handle one-to-many mappings (e.g., one GO term can be assigned to multiple GO Slim terms).

```bash
woltka tools collapse \
  --input output/truth.uniref.tsv \
  --map function/go/goslim.tsv.xz \
  --names function/go/name.txt.xz \
  --output truth.goslim.tsv
```
