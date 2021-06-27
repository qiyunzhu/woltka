# Output files

The main output file(s) of Woltka are **feature tables** (a.k.a., **profiles**), in which column headers are sample IDs and row headers are feature IDs, and each cell value represents the _count_ of a particular feature in a particular sample. For example:

Feature ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 | Name | Rank
--- | --- | --- | --- | --- | --- | ---
561 | 20 | 45 | 5 | 15 | Escherichia | genus
570 | 135 | 110 | 75 | 150 | Klebsiella | genus
590 | 0 | 10 | 525 | 5 | Salmonella | genus
543 | 5 | 10 | 0 | 12 | Enterobacteriaceae | family
... |

## Contents

- [Output filepath](#output-filepath)
- [File formats](#file-formats)
- [Output read maps](#output-read-maps)
- [Table utilities](#table-utilities)

## Output filepath

Parameter `--output` or `-o` is to tell Woltka where to save the output feature table(s).

The number of feature tables generated per run is dependent on the number of ranks (supplied with `--rank` or `-r`, see [classification](hierarchy.md) for details). If there is only one rank, the value of the parameter is the actual path of the output file. For example:

```bash
woltka classify ... -r genus -o output.biom
```

If there are multiple ranks, this parameter points to a directory, where multiple feature tables with the filename pattern `<rank>.biom` will be stored. For example:

```bash
woltka classify ... -r phylum,genus,species -o .
```

## File formats

### BIOM

The default output feature table format is [**BIOM**](http://biom-format.org/), a general format for efficiently storing and accessing omics data tables. A feature table file will have the extension name `.biom`.

A BIOM table not only stores feature counts per sample, but also allows appending feature metadata to it. Woltka optionally appends the following three metadata columns to a BIOM table:

- **Name**: The name of each classification unit, as defined in the file provided by `--names` (see [classification](hierarchy.md) for details).

  Alternatively, instead of appending the name column, one may choose to replace unit IDs with their names, by using the `--name-as-id` flag. IDs without names will remain IDs. Note however, this will throw an error if multiple IDs correspond to the same name!

- **Rank**: The rank of each classification unit, if available. This column is disabled by default, and it can be enabled by using the `--add-rank` flag.

- **Lineage**: A string of units ordered from high to low and joined by semicolon. This is also known as the "[Greengenes](https://greengenes.secondgenome.com/)"-style lineage string, which is widely used in bioinformatics. For example:

  ```Viruses;Microviridae;Bullavirinae;Sinsheimervirus;Escherichia_virus_phiX174```

  This column is disabled by default, and it can be enabled by using the `--add-lineage` flag.

  The aforementioned `--name-as-id` flag also affects how the lineage strings are composed (IDs or names).

A BIOM table can be directly parsed by many bioinformatics programs. It is typically not necessary to convert it back to plain text format (i.e., a tab-separated file). However if you do want to do so, the following command with do the job:

```bash
biom convert --to-tsv -i input.biom -o output.tsv
```

Optionally, one can append a metadata column to the right of the table. For example:

```bash
biom convert --to-tsv -i input.biom -o output.tsv --header-key Name --tsv-metadata-formatter naive
```

Additionally, one can export all metadata as a separate table, with:

```bash
biom export-metadata -i input.biom --observation-metadata-fp metadata.tsv
```

### TSV

In addition to BIOM, Woltka supports directly outputing the old-style tab-separated files (**TSV**). This can be toggled as follows:

If the output filepath (`--output`) is a single file (i.e., `--rank` has only one value), simply name it as `xxx.biom` to write in BIOM format, or anything else (e.g., `xxx.tsv`) to write in TSV format.

Alternatively, use `--to-biom` or `--to-tsv` to force a specific output format. This is useful when there are multiple ranks and the output filepath is a directory. Example:

```bash
woltka classify ... -r phylum,genus,species -o . --to-tsv
```

Metadata (see above) will be appended as extra columns to the right of the table.

Likewise, one can convert a TSV file (without metadata) to BIOM format:

```bash
biom convert --to-hdf5 -i input.tsv -o output.biom
```


## Output read maps

During the processing of each alignment file, Woltka can optionally write to disk a map of reads to classification units. These "read maps" are informative when the user wants to know how exactly individual reads are assigned. They are useful in a series of downstream applications. For example, Woltka's stratification function can take read maps generated using one classification system (e.g., taxonomy), and stratify assignments based on another (e.g., function). See [stratification](stratify.md) for details.

The format of a read map is usually as simple as "read \<tab\> classification unit". If ambiguous assignment is allowed (`--ambig`, see [classification](hierarchy.md) for details), one or multiple classification units will be joined by comma, and counted by a number following colon if the same classification unit occurs multiple times. Example:

```
R01 Ecoli
R02 Ecoli,Strep
R03 Ecoli:2,Kleb:3,Cdiff
...
```

The flag `--name-as-id` also applies to the read maps (see [above](#file-formats)).

To enable writing of read maps, add `--outmap <directory>` to the command line. Woltka will write one read map file per sample, named after the sample ID, to the directory. If there are multiple ranks (`--rank`) to which reads are classified, Woltka will write read maps at each rank in a separate subdirectory named by the particular rank.

Read maps are typically large -- comparable to the original alignment files, since each read occupies one line. Therefore, Woltka by default compresses them using the gzip algorithm (extension: `.gz`). One can select compression method using the `--outmap-zip` parameters. Choices are `none`, `gz`, `bz2` and `xz`.

For example, with command:

```bash
woltka classify ... -s S1,S2,S3 -r genus --outmap maps
```

In directory `maps`, there will be three read map files: `S1.txt.gz`, `S2.txt.gz` and `S3.txt.gz`.

With command:

```bash
woltka classify \
  ...
  --samples S1,S2,S3 \
  --rank phylum,genus,species \
  --outmap outmap_dir \
  --outmap-zip xz
```

In `outmap_dir`, there will be three subdirectories: `phylum`, `genus` and `species`, eaching holding three read map files: `S1.txt.xz`, `S2.txt.xz` and `S3.txt.xz`.


## Table utilities

Woltka provides several utilities under the `tools` menu for table manipulation (both BIOM and TSV are supported and automatically recognized). They are:

- [**Collapse**](collapse.md): Collapse a profile based on a source-to-target feature mapping; supporting many-to-many relationships.
- [**Coverage**](coverage.md): Calculate per-sample coverage of feature groups, such as the completeness of metabolic pathways or core gene sets.
- [**Normalize**](normalize.md): Normalize a profile based on feature sizes or into relative abundances.
- [**Filter**](filter.md): Filter features by per-sample abundance.
- [**Merge**](merge.md): Merge multiple profiles into one.
