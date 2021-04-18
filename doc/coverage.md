# Feature group coverage

The **coverage** command calculates the coverage -- percentage of features present in each sample over a pre-defined group of features -- of a profile.

```bash
woltka tools coverage -i input.biom -m mapping.txt -o output.biom
```

A typical use case is to assess the likelihoods of presence of **metabolic pathways** in each organism or community. Because a pathway consists of _multiple_ chemical **reactions** or functional **genes** connected to each other, the presence of some of them (even with high abundance) in the sample does not necessarily suggest that the entire pathway is viable. Only when all or a large proportion of them are found can we be more confident about this hypothesis.

In this example, the input profile ([sample](../woltka/tests/data/output/truth.metacyc.tsv)) is a table of **genes**:

Feature ID | Sample 1 | Sample 2 | Sample 3 | Sample 4
--- | --- | --- | --- | ---
_plsC_ | 51 | 49 | 113 | 34
_fruK_ | 83 | 128 | 160 | 41
_panE_ | 0 | 53 | 0 | 39
_leuA_ | 111 | 262 | 232 | 77
... |

The mapping file ([sample](../woltka/tests/data/function/metacyc/pathway_mbrs.txt)) defines the member features (**genes**) of each feature group (**pathway**) (each line can have arbitrary number of fields; field delimiter is \<tab\>):

| | | | | | | |
|-|-|-|-|-|-|-|
| Asparagine biosynthesis | _asnB_ | _aspC_ |
| Biotin synthesis | _bioA_ | _bioB_ | _bioD_ | _bioF_ |
| NAD biosynthesis II | _hel_ | _nudC_ | _nadN_ | _pnuE_ | _nadR_ | _nadM_ |
| pyruvate decarboxylation | _aceE_ | _aceF_ | _lpd_ |
| ... |

The output file ([sample](../woltka/tests/data/output/truth.metacyc.coverage.tsv)) is a table of coverage values (percentages) per sample per feature group (**pathway**):

Feature ID | Sample 1 | Sample 2 | Sample 3 | Sample 4
--- | --- | --- | --- | ---
Biotin synthesis | 50.0 | 50.0 | 25.0 | 37.5
GDP-D-rhamnose biosynthesis | 20.0 | 80.0 | 20.0 | 80.0
L-glutamine degradation I | 100.0 | 100.0 | 50.0 | 0.0
Sucrose biosynthesis I | 20.0 | 20.0 | 20.0 | 20.0
... |

**Note**: The "coverage" computed by Woltka is distinct from those by HUMAnN2 (whether the pathway is present) and HUMAnN3 (how likely the pathway is present).

## Parameters

### Presence / absence

With parameter `--threshold` or `-t` followed by a percentage (e.g., `80`), the output coverage table will display binary results, with "**1**" representing coverage above or equal to this threshold and "**0**" being coverage below this threshold.

### Feature count

With flag `--count` or `-c`, the program will report the number of member features of a group present in a sample, instead of the percentage. Note: This will override `--threshold`.

### Feature group names

One can supply a mapping of feature groups to their names by `--names` or `-n`, and these names will be appended to the coverage table as a metadata column ("Name").


## Considerations

The coverage command will treat any feature count -- as low as **1** -- as the evidence of the feature's presence. False positives may be introduced if the profile has many noises. One may consider **filtering** the profile prior to running this command. Woltka provides a per-sample feature abundance [filtering](filter.md) function, in addition to the multiple filtering functions implemented in the QIIME 2 plugin [feature-table](https://docs.qiime2.org/2020.11/plugins/available/feature-table/).
