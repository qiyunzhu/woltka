# Data for Unit Test

Test datasets, including reference genome and taxonomy databases, query sequencing data are based on a a collection of 107 NCBI-defined "reference" bacterial genomes. The database is available for download at:

## Align

**Sequence alignment** results. Here "alignment" refers to the operation of aligning short DNA sequences ("reads") against reference genome sequences. The query sequences are five samples (S01 to S05) of 150 bp paired-end reads simulated using [CAMISIM](https://github.com/CAMI-challenge/CAMISIM/). Three aligners were used: [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Centrifuge](https://ccb.jhu.edu/software/centrifuge/), and [BURST](https://github.com/knights-lab/BURST). In particular, `bowtie2` used the Bowtie2 program default, whereas `bt2sho` used the Bowtie2 parameter setting recommended by SHOGUN. The mappings between reads and genomes are provided by each alignment file. In addition, `truth` is the ground-truth mapping derived from CAMISIM output.

## Taxonomy

**Reference taxonomy** of the 107 bacterial genomes. The taxonomy is provided in multiple formats:

- `g2lineage.txt`: Greengenes-style lineage strings.
- `taxdump`: NCBI taxdump-style database files (`nodes.dmp` and `names.dmp`).
- `g2tid.txt`: Genome ID to NCBI taxonomy ID mapping file.
- `rank_names.tsv` and `rank_tids.tsv`: Taxon name and TaxID at each of the seven standard taxonomic ranks for each genome.

In addition, `nucl` contains the mappings from nucleotide sequence accessions (instead of their host genomes) to taxonomy.

## Phylogeny

- `tree.nwk`: maximum likelihood phylogenetic tree of the 107 bacterial genomes built using RAxML and midpoint-rooted.
