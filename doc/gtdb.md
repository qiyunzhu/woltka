# Working with GTDB

[**GTDB**](https://gtdb.ecogenomic.org/) (Genome Taxonomy Database) ([Parks et al., 2018](https://www.nature.com/articles/nbt.4229)) is a standardized taxonomy system for bacteria and archaea. It is based on the phylogenetic trees of genomes, hence reflecting the evolutionary relationships among microorganisms better, compared to conventional taxonomy. Our research shows that GTDB shares high consistency with the [WoL](https://biocore.github.io/wol/) phylogeny (which was built using more marker genes and more expensive methods) ([Zhu et al., 2019](https://www.nature.com/articles/s41467-019-13443-4)).

This tutorial will be based on GTDB [**R04-RS89**](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/) (version 4, indexed to RefSeq release 89), which was the current release at the time of writing. This release contains 143,512 bacterial and 2,392 archaeal organisms in its taxonomy, or 23,458 bacterial and 1,248 archaeal genomes in its collection of actual genome sequences. These 24,706 "**species clusters**" cover the complete GTDB taxonomic framework ([Parks et al., 2020](https://www.nature.com/articles/s41587-020-0501-8)).

This tutorial will explain various strategies of utilizing GTDB taxonomy / phylogeny in a meta'omics analysis with Woltka.


## Contents

- [Woltka directly parses GTDB taxonomy](#woltka-directly-parses-gtdb-taxonomy)
- [Reformat GTDB as taxdump style](#reformat-gtdb-as-taxdump-style)
- [Build GTDB genome database](#build-gtdb-genome-database)
- [Use GTDB phylogeny](#use-gtdb-phylogeny)
- [Use WoL database with GTDB taxonomy](#use-wol-database-with-gtdb-taxonomy)


## Woltka directly parses GTDB taxonomy

The GTDB taxonomy ([`bac120_taxonomy_r89.tsv`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_taxonomy_r89.tsv) and [`ar122_taxonomy_r89.tsv`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv)) is provided in the format of **Greengenes-style lineage strings**. For example: _Escherichia coli_ O157:H7 str. Sakai (NCBI RefSeq assembly: [GCF_000008865.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000008865.1/)) is classified as:

```
RS_GCF_000008865.1	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia dysenteriae
```

This format can be directly parsed by Woltka:

```bash
woltka classify --lineage bac120_taxonomy_r89.tsv --lineage ar122_taxonomy_r89.tsv ...
```

Note: The [GTDB-Tk](https://github.com/Ecogenomics/GtdbTk) database ([`gtdbtk_data.tar.gz`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz)) already provides a combined Bacteria and Archaea taxonomy file: `gtdb_taxonomy.tsv`.


## Reformat GTDB as taxdump style

Alternatively, you may consider converting GTDB taxonomy into the format of NCBI taxdump. This will allow you to connect GTDB to more bioinformatics tools such as BLAST and Kraken2 which favor the taxdump format. This will also make the Woltka output slimmer, if that's your preference. For example, instead of the long lineage strings like `d__Bacteria;p__Proteobacteria;...`, you will get taxon names like `Escherichia` as the feature name.

Three good characteristics of the GTDB taxonomy (in addition to many others) are that:

1) All lineages have exactly seven ranks (`d`omain, `p`hylum, `c`lass, `o`rder, `f`amily, `g`enus and `s`pecies).
2) All ranks are filled (i.e., no "gaps" such as `g__;` as multiple other systems do).
3) At each rank, all taxon names are unique (but duplicate names are allowed at different ranks).

These characteristics make it safe to collapse GTDB lineages into taxonomic units and index them in a way that resembles NCBI taxdump. (Meanwhile, in case you are thinking about whether the method explained here can also be applied to other taxonomies, you need to keep these in mind.)

We provide a script: [`gtdb_to_taxdump.py`](https://biocore.github.io/wol/code/scripts/gtdb_to_taxdump.py) to automate this conversion:

```bash
gtdb_to_taxdump.py bac120_taxonomy_r89.tsv ar122_taxonomy_r89.tsv
```

If you start with `gtdb_taxonomy.tsv` from [GTDB-Tk](https://github.com/Ecogenomics/GtdbTk), do this instead:

```bash
gtdb_to_taxdump.py gtdb_taxonomy.tsv
```

This command will generate three files: `nodes.dmp`, `names.dmp` and `taxid.map`. They constitute an NCBI-style taxonomy database (to learn more please refer to the official [NCBI taxonomy database](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)), in which every taxon is given a taxonomy ID. This ID is an incremental integer, starting from "1" at the root, and increasing from higher ranks to lower ranks. So domains Bacteria and Archaea receive "2" and "3", respectively, followed by phyla, so on so forth.

With these output files, you can simply do:

```bash
woltka classify --map taxid.map --nodes nodes.dmp --names names.dmp ...
```


## Build GTDB genome database

Now we know that Woltka can handle GTDB taxonomy, but what about the input alignment files? Are there some considerations when we build the GTDB genome database for specific aligners (instead of directly running [GTDB-Tk](https://github.com/Ecogenomics/GtdbTk)) so that the resulting alignment files can be fed into Woltka, together with the GTDB taxonomy?

As mentioned above, the GTDB data release provides nucleotide sequences of 24,706 bacterial and archaeal genomes. These sequences are in archive [`gtdb_r89_rep_genomes.tar.gz`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdb_r89_rep_genomes.tar.gz). Unzip and move the files under `bacteria/` and `archaea/` into one directory. Now you have 24,706 gzipped Fasta files, with the filename pattern as:

- NCBI RefSeq genomes: `RS_GCF_000123456.1_genomic.fna.gz`
- NCBI GenBank genomes: `GB_GCA_000654321.1_genomic.fna.gz`
- MAGs found by the GTDB team: `UBA12345_genomic.fna.gz`

We provide a script: [`genomes_for_db.py`](https://biocore.github.io/wol/code/scripts/genomes_for_db.py) to convert multiple genome sequence files into a single Fasta file.

```bash
genomes_for_db.py --input input_dir --output output.fna --ext _genomic.fna.gz --concat --gap N*20
```

This command will removing intra-sequence line breaks, concatenate multiple sequences in the same file (e.g., chromosomes or scaffolds of the same genome) with a linker of 20 "N" characters, and merge all genomes into one multi-Fasta file, in which sequence IDs are genome IDs.

The output file `output.fna` can then be used for database building with choice of tools. The taxonomy files discussed above will match the results of sequence alignment using this database. For examples:

Bowtie2:

```bash
bowtie2-build --seed 42 --threads 32 output.fna GTDBr89
```

BLASTn (with taxonomy integrated, see [above](#reformat-gtdb-as-taxdump-style)):

```bash
makeblastdb -in output.fna -out GTDBr89 -dbtype nucl -parse_seqids -taxid_map taxid.map
```

Note: The same set of genomes can also be found in the [GTDB-Tk](https://github.com/Ecogenomics/GtdbTk) database, under `fastani/database/`. However these filenames do not have the prefix `RS_` or `GB_`. Therefore, if you choose to use these files, you will need to modify the taxonomy table to match genome IDs:

```bash
cat taxonomy/gtdb_taxonomy.tsv | sed 's/^RS_\|GB_//' > gtdb_taxonomy_rev.tsv
```


## Use GTDB phylogeny

GTDB is a phylogeny-based taxonomy system, which means that instead of staying with the derived, hard-coded taxonomic ranks, you may directly exploit the original phylogenetic trees. However there are some tricks to be considered:

The "GTDB phylogeny" actually consists of two trees: one for Bacteria ([`bac120_r89.sp_labels.tree`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_r89.sp_labels.tree)) and the other for Archaea ([`ar122_r89.sp_labels.tree`](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_r89.sp_labels.tree)).

```bash
woltka classify --newick bac120_r89.sp_labels.tree --newick ar122_r89.sp_labels.tree ...
```

Woltka will automatically join the two trees into one at their roots. This is operational, but note that the two GTDB trees are NOT rooted based on evolutionary origins, so it is your decision whether this merging is fair.


## Use WoL database with GTDB taxonomy

If you decide to stick to the WoL genome database (we appreciate that!) but also want to embrace the GTDB taxonomy, we have a ready-to-go solution for you:

We provide original and curated NCBI and GTDB [taxonomies](https://biocore.github.io/wol/data/taxonomy/) based on the WoL phylogeny. For GTDB specifically, you may download this file: [`lineages.txt.bz2`](https://github.com/biocore/wol/raw/master/data/taxonomy/gtdb/curation/lineages.txt.bz2), which contains curated GTDB lineage strings for WoL genomes.

Similar to above, you can simply do:

```bash
woltka classify --lineage lineages.txt.bz2 --rank phylum,genus,species ...
```

As you would expect, the resulting taxonomic profiles will speak GTDB!
