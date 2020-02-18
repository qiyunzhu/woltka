# q2-woltka

**q2-woltka** is a QIIME 2 plugin for Woltka, the metagenome data analysis tool.


## Installation

Requires: QIIME 2 2019.1 or above.

First, install the standalone Woltka program:

```bash
pip install git+https://github.com/qiyunzhu/woltka.git@dev
```

Second, update QIIME 2 cache and the plugin will be automatically loaded into QIIME 2.

```bash
qiime dev refresh-cache
```

Then one can launch the program by executing:

```bash
qiime woltka
```

## Example usage

Woltka ships with small test datasets under this directory:

```
<program_dir>/woltka/tests/data
```

One can execute the following commands to make sure that Woltka functions correctly, and to get an impression of the basic usage of Woltka.

gOTU table generation:

```bash
qiime woltka gotu --p-align-dir align/bowtie2 --o-table table.qza
```

The output file `table.qza` contains a BIOM table (QIIME artifact type: `FeatureTable[Frequency]`) in which columns are samples and rows are gOTUs. It can then be analyzed using the classical QIIME 2 workflow.
