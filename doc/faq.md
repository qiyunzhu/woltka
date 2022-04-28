# Frequently asked questions


## Computing

### Is Woltka deterministic or stochastic?

Woltka is **deterministic**. Given the same input files and parameters, it always produces the identical output files. There is no "seed" parameter.

### Is Woltka exhaustive/exact/optimal? Or is it approximate/heuristic?

The former. Woltka **exhaustively** captures all valid matches from the alignment file(s).

### Are Woltka results consistent across versions?

To date, all Woltka versions (0.1.0 to 0.1.4) generate **identical** output files given the same setting. Later versions are more efficient and have more features, though.

### How many CPU cores does Woltka use?

Woltka works the best with two CPU cores: one for file decompression and the other for classification. This happens automatically. See [here](perform.md#keep-external-decompressors-on) for details.

### Does Woltka support multiprocessing?

Not at the moment. But you can run multiple Woltka instances on different subsets of samples and merge results. See [here](perform.md#run-separate-jobs-and-merge-results) for details.


## Input files

### Does Woltka support gzipped alignment files?

Yes. All input files for Woltka (alignments and databases) can be supplied as compressed in gzip, bzip2 or xz formats. Woltka will automatically recognize and process them.

### Does Woltka support [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) and [CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format)) formats?

Not out-of-the-box. But you can use SAMtools to extract BAM/CRAM files and directly "pipe" into Woltka, like this ("-" represents stdin):

```bash
samtools view input.bam | woltka classify -i - -o output.biom
```

### Does Woltka support [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) format?

Not out-of-the-box. But you can use the following AWK trick to convert a PAF file into mock BLAST format and feed into Woltka. There will be no percent identity, e-value or bit score, but Woltka doesn't need them anyway.

```bash
cat input.paf | awk -v OFS="\t" '{print $1,$6,0,$11,0,0,$3+1,$4,$8+1,$9,0,$12}' | woltka classify -i - -o output.biom
```

### I ran `woltka classify -i input.fastq ...`, and got an error saying it cannot determine alignment file format. Why?

Woltka takes alignment files as input, NOT original sequencing data (FASTQ, FASTA, etc.). You need to perform alignment on the sequencing data by yourself, such as:

```bash
bowtie2 -x db -f input.fastq -S output.sam
```

Then feed the resulting alignment(s) into Woltka.

```bash
woltka classify -i output.sam ... -o output.tsv
```

See [here](align.md) for details.

### I ran `woltka classify -i S01.sam ...`. The output feature table has a single column with an empty header. Why?

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

### Can Woltka report cell values in the units of CPM, RPK, and TPM?

Yes. See [here](normalize.md) for methods.

### I ran Woltka separately on multiple subsets of data. Can I merge the results?

Yes. The `woltka tools merge` command is for you. See [here](merge.md) for details.

### Can Woltka report taxon names instead of TaxIDs when using NCBI taxonomy?

Yes. Add `--name-as-id` to the command? See [here](output.md) for details.


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
