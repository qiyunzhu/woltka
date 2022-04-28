# Working with RefSeq

The NCBI Reference Sequence Database (**RefSeq**) (https://www.ncbi.nlm.nih.gov/refseq/) ([O'Leary et al., 2016](https://academic.oup.com/nar/article/44/D1/D733/2502674)) is a database of curated, non-redundant biological sequences. Genome sequences are released at the [RefSeq genomic FTP](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) server. It is among the largest public sources of genomic data.

We provide a Python script: [**refseq_build.py**](https://github.com/qiyunzhu/utils/blob/main/refseq_build.py), which can automatically sample latest genomic data hosted at the RefSeq database, download them, and reformat them such that Woltka (or other programs) can utilize the resource. The script has various options for customizing the database to address the best need of the user.


## Contents

- [Protocols](#protocols)
- [Customization](#customization)
- [Database files](#database-files)
- [Using the database](#using-the-database)
- [Adding extras](#adding-extras)

## Protocols

RefSeq is large, therefore one typically cannot exhaustively include all RefSeq genomes in one database. Sampling is necessary.

As of 2022-01-01, RefSeq contains 248,391 genomes. This number becomes 1,170,106 if GenBank genomes are also considered (excluding redundant ones). The numbers of genomes sampled using individual protocols are listed below. These numbers will likely grow in the future, considering [how fast RefSeq is growing](https://www.ncbi.nlm.nih.gov/refseq/statistics/) every day.

The script [**refseq_build.py**](https://github.com/qiyunzhu/utils/blob/main/refseq_build.py) by default samples genomes from five categories: **archaea**, **bacteria**, **fungi**, **protozoa**, and **viral**. They are most relevant to microbiome studies. But the user can customize this using the `--cats` parameter.

### Basic (3,522 genomes)

```bash
refseq_build.py -p basic -o .
```

The basic protocol samples one **complete** genome per **genus** (defined by NCBI). It is equivalent to the following command.

```bash
refseq_build.py --complete --sample 1 --rank genus -o .
```

Suitable for coarse classification. Also suitable for practice and demonstration purpose.

### Standard (15,403 genomes)

```bash
refseq_build.py -p standard -o .
```

Equivalent to:

```bash
refseq_build.py --sample 0 --reference --represent -o .
```

The standard protocol simply collects all NCBI-defined **reference and representative** genomes (a.k.a., [RefSeq Select](https://www.ncbi.nlm.nih.gov/refseq/refseq_select/)). They are useful for exploring inter-species diversity of microbiomes found in common habitats.

We **recommend** this protocol for general purpose microbiome classification.

### Extended (29,603 genomes)

```bash
refseq_build.py -p extended -o .
```

Equivalent to:

```bash
refseq_build.py --sample 1 --rank species_latin --above --reference --represent --typemater -o .
```

The extended protocol samples one genome per species that has a Latinate name, as well as higher ranks (genus to phylum, in case something is missing). It also includes all reference, representative, and [type material](https://academic.oup.com/nar/article/43/D1/D1086/2438106) genomes.

- The sampling priority is that if a species has a genome that belongs to the three special categories, then select it; if not, then sort the genomes by the order of complete genome > scaffolds > contigs, and select the first one.

### Extreme (95,244 genomes)

```bash
refseq_build.py -p extreme -o .
```

Equivalent to:

```bash
refseq_build.py --genbank --sample 1 --rank species --above --reference --represent --typemater -o .
```

With a reference genome database of this size, most current sequence alignment tools will struggle, however alignment-free methods (such as [Mash](https://github.com/marbl/Mash)) may still work. The user should decide whether it is worth.


Finally, just for reference, if one wants to pull down the entire NCBI, the command is:

```bash
refseq_build.py --cats all --block "" --genbank -o .
```


## Customization

The **refseq_build.py** has high customizability, with which one may sample genomes that are most relevant to specific research goals while retaining computational feasibility (that is, not too many genomes, just relevant ones).

### More examples

```bash
refseq_build.py -c bacteria,archaea -o . -t 1117,766 -e
```

This will include bacterial and archaeal genomes, but exclude ones classified as Cyanobacteria (TaxID: [1117](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1117)) and Rickettsiales (TaxID: [766](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=780)). Useful for ruling out susceptible chloroplastic or mitochondrial signals.

```bash
refseq_build.py --genbank --sample 10 --rank phylum -o .
```

This will include GenBank genomes, and select up to 10 genomes per phylum. This small database will maximize the inclusion of the 100+ mostly uncultivated, understudied candidate phyla which has less curated genomes. Useful for exploring deep branches of life.

```bash
refseq_build.py -g gids.txt -o .
```

This will only download genomes specified in the file `gids.txt`. Useful for controlled tests.

### Break and resume

Should any of the download steps be interrupted by e.g., a network failure, one can resume the downloading process by re-executing the same command. The program will skip the already downloaded files in this new run. In some instances, one may need to manually remove the last file from the failed run (because that file may be corrupt), before re-running the program.

If one wants to overwrite downloaded files (e.g., upgrading), add `--overwrite` to the command.

### Manual downloading

One may want to download genomes manually in a more controled manner, instead of letting the script running for hours to days to retrieve them one after another before moving to the next step. In this case, add `--manual` to the command, and the program will generate `urls.txt`, a list of URLs of the sampled genomes, and quit.

Then one can choose the most appropriate method to download them. For example, one may use the [rsync protocol](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#protocols), as recommended by NCBI:

```bash
while read url; do
  for x in fna gff; do
    rsync -Ltq rsync${url#https}/${url##*/}_genomic.$x.gz download/genomic_${x}/
  done
done < urls.txt
```

After all genomes are downloaded, one may restart the program without `--manual`, and the program will take the downloaded files and move to the next step.


## Database files

After running **refseq_build.py**, the output directory has these files and subdirectories:

Name | Description
--- | ---
`metadata.tsv` | Metadata of sampled genomes.
`db.fna` | A multi-Fasta file containing all genome sequences. It is ready to be [indexed](align.md#alignment-with-bowtie2) by specific alignment tools.
`taxid.map` | Genome ID to taxonomy ID mapping file.
`nodes.dmp` and `names.dmp` | Shrinked NCBI taxonomy database that only contain taxonomic information of sampled genomes and higher-rank units.
`lineages.txt` | Genome ID to Greengenes-style lineage strings. Only the seven standard ranks are preserved. An alternative to the NCBI style.
`length.map` | Genome ID to total sequence length (bp) mapping. Useful for [normalization](normalize.md#normalization-by-subject-size).
`coords.txt` | Coordinates of genes (CDS sequences, as identified by their locus tags) on their host genomes. Useful for [coord-matching](ordinal.md).
`gene.map`, `protein.map` and `product.map` | Gene annotations, including gene symbols, protein IDs and product descriptions, respectively. Useful for functional profiling.
`download/` | Directory containing raw files downloaded from the NCBI server. This includes the assembly summary of all RefSeq genomes (`assembly_summary_refseq.txt`), the full-size taxonomy database (`taxdump`), the genome ID lists of individual categories (`cats/`), the genomic DNA sequences (`genomic_fna/`), and the genome annotations (`genomic_gff`).

You may compress any of these files to save disk space. Woltka can [handle compressed files](perform.md#compress-alignment-files). You may delete the `download/` directory if you no longer need the raw files.


## Using the database

The procedures discussed in other documents also apply here. For examples:

Index the genome sequences using Bowtie2:

```bash
bowtie2-build --threads 8 db.fna db
```

Perform sequence alignment on individual samples:

```bash
while read $id; do 
  bowtie2 --very-sensitive -p 8 -x db -1 $id.R1.fq -2 $id.R2.fq -S bt2out/$id.sam
done < sample.list
```

Generate an OGU table:

```bash
woltka classify -i bt2out -o output.biom
```

Taxonomic profiling by species:

```bash
woltka classify -i bt2out -m taxid.map --nodes nodes.dmp --names names.dmp -r species -o output.biom
```

Matching reads with genes:

```bash
woltka classify -i bt2out -c coords.txt -o output.biom
```

Functional profiling by gene symbol: 

```bash
woltka classify -i bt2out -c coords.txt -m gene.map -o output.biom
```


## Adding extras

RefSeq is probably the largest database, but it isn't the most feature-rich one. Some components that are useful in downstream analyses are not readily available. We would like to share our thoughts on what you can add to the resource.


### Phylogenetic tree

Unlike our [WoL database](wol.md), RefSeq doesn't ship with a reference phylogeny. Therefore one cannot perform [phylogeny-aware analyses](https://journals.asm.org/doi/10.1128/msystems.00167-22) out-of-the-box.

To build a phylogenetic tree, the tough (but high-quality) way is to follow the [WoL protocol](https://biocore.github.io/wol/protocols/). However, there is a convenient alternative: using [PhyloPhlAn](https://github.com/biobakery/phylophlan) ([Asnicar et al., 2020](https://www.nature.com/articles/s41467-020-16366-7)) (which some of us also contributed to). With PhyloPhlAn, building a phylogenomic "tree of life" is as simple as:

```bash
phylophlan -i genomic_fna/ -d phylophlan -f tol.cfg --diversity high --fast -o output_tol --nproc 8
```

It should be noted that this method does not handle eukaryotes (e.g., protists and fungi). Meanwhile, it is not easy (if ever possible) to include viruses in the tree of life. Therefore, to build a tree like this, you will have to skip these organism groups.

### Function catalogs

The RefSeq genome annotation files (.gff) do not provide direct mapping from genes to higher functional catalogs. However, one may obtain these mappings by exploring external database resources and data mining technologies.

For example, locus "M6_RS00005" has a gene symbol "_dnaA_" and a RefSeq protein ID "[WP_012657571.1](https://www.ncbi.nlm.nih.gov/protein/WP_012657571.1/)". The RefSeq protein record contains mappings to Gene Ontology (GO) entries (such as [GO:0005524 - ATP binding](http://amigo.geneontology.org/amigo/term/GO:0005524)). These links may be batch-retrieved using [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/).

For another example, the [UniProt](https://www.uniprot.org/downloads) data release contains "ID mapping". This includes mapping from UniProt entries to RefSeq protein IDs. In the example above, WP_012657571.1 can be mapped to [B9DSN7](https://www.uniprot.org/uniprot/B9DSN7). With this file, one can build a reverse mapping to convert RefSeq proteins to UniProt entries. Then one can use other "ID mappings" in the UniProt data release to navigate to other functional catalogs.
