# Frequently asked questions


## Computing

### Is Woltka deterministic or stochastic?

Woltka is **deterministic**. Given the same input files and parameters, it always produces the identical output files. There is no "seed" parameter.

### Is Woltka exhaustive/exact/optimal? Or is it approximate/heuristic?

The former. Woltka **exhaustively** captures all valid matches from the alignment file(s).

### How many CPU cores does Woltka use?

Woltka works the best with two CPU cores (threads): one for decompression and the other for classification. This is automatically determined.


## Input files

### Can Woltka parse compressed files?

Yes. All input files for Woltka (alignments and databases) can be supplied as compressed in gzip, bzip2 or xz formats. Woltka will automatically recognize and parse them.

### I ran `woltka classify -i input.fastq -o output.tsv`, and got an error saying it cannot determine alignment file format. Why?

Woltka takes alignment files as input, NOT original sequencing data (FASTQ, FASTA, etc.). You need to perform alignment on the sequencing data by yourself, such as:

```bash
bowtie2 --very-sensitive -x db -f input.fastq -S output.sam
```

Then feed the resulting alignment(s) into Woltka.

```bash
woltka classify -i output.sam ... -o output.tsv
```

### I ran `woltka classify -i S01.sam ... -o S01.tsv`. The output feature table has a single column with an empty header. Why?

Woltka is designed to deal with multiple samples at once. If the input is a single alignment file, Woltka automatically treats it as a multiplexed file and attempts to extract individual samples out of it. But in case it is actually not -- Woltka will leave the sample ID empty.

If you really wants Woltka to process one sample at a time, you can do this (aside from manually adding sample ID to the output table):

```bash
woltka classify --no-demux -i S01.sam ... -o S01.tsv
```

The `--no-demux` flag will tell Woltka not to try to demultiplex the alignment file. Instead, it will use the filename `S01` as the only sample ID.

See [here](input.md#demultiplexing) for details.


## Output files

### Are cell values in the output feature table absolute or relative?

By default, cell values are feature _frequencies_, i.e., the numbers of sequencing reads assigned to individual features. Therefore they are **absolute** abundance.

Woltka provides multiple normalization features. For example, if you want to get **relative** abundance (fraction) instead of absolute abundance, you can do:

```bash
woltka classify ... --frac
```

See [here](normalize.md) for details.

### I ran Woltka separately on multiple subsets of data. Can I merge the results?

Yes. The `woltka tools merge` command is for you. See [here](merge.md) for details.


## Hierarchies

### What if a read is aligned to multiple reference sequences?

By default, Woltka will keep all matches and divide them by the number of matches. Meanwhile, Woltka lets the user choose from multiple alternative behaviors.

See [here](classify.md#ambiguous-assignment) for details.

### What if some subjects are not found in the classification system?

They will not be counted in the output feature table, unless you specify the `--unassigned` flag, in which case an extra feature `unassigned` will be appear.

See [here](classify.md#unassigned-sequences) for details.


## Stratification

### I attempted to collapse a stratified feature table, and the output is empty?

In this case, you will need to specify which field of a stratified feature should be collapsed, using parameter `--field` followed by field index (starting from 1), otherwise the program will try to collapse the entire feature.

See [here](collapse.md#stratification) for details.
