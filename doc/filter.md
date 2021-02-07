# Per-sample filtering

The **filter** command filters each feature in each sample based on the absolute or relative abundance of that feature in that particular sample. For example, the following command will drop features that are less than 0.01% abundant in each sample:

```bash
woltka tools filter -i input.biom -o output.biom --min-percent 0.01
```

This function is especially useful in shotgun metagenomics, where very-low-abundance false positive assignments are prevalent and causing biases in downstream analyses ([Ye et al, 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5)).
