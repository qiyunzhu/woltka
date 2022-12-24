# Per-sample filtering

The **filter** command filters each feature in each sample based on the absolute or relative abundance of that feature in that particular sample. For example, the following command will drop features that are less than 0.01% abundant in each sample:

```bash
woltka filter -i input.biom -o output.biom --min-percent 0.01
```

This function is especially useful in shotgun metagenomics, where **very-low-abundance false positive assignments** are prevalent and can cause biases in downstream analyses ([Ye et al, 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5)).

Whether one should filter the feature table depends on the quality of input alignments and the methods to be used in the downstream analysis. As a rule of thumb, methods that consider 1) abundance of individual features and 2) relations among features (such as a phylogenetic tree) are more robust against low-abundance false-positives. One example is **weighted UniFrac**. If the downstream analysis relies on these methods, filtering of the feature table is less necessary.
