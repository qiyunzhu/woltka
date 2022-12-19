# Working with KEGG

**KEGG** (https://www.genome.jp/kegg/) ([Kanehisa et al., 2021](https://academic.oup.com/nar/article/49/D1/D545/5943834)) is a classical database of biological functions. It provides a well-organized hierarchical identification system, such as orthologies (K), modules (M), reactions (R), compounds (C), pathways, diseases and more.

Whereas the FTP access to KEGG is limited to users with a subscription, the mapping of UniRef entries to KEGG orthology (KO) entries is freely available from the [UniProt](https://www.uniprot.org/downloads) data release. From this point on, we provide a Python script: [**kegg_query.py**](https://github.com/qiyunzhu/utils/blob/main/kegg_query.py) to automatically retrieve higher-level classification information of a given KO list or table from the KEGG server. This is made possible using the official [KEGG REST API](https://www.kegg.jp/kegg/rest/), which is freely available to academic users (however restrictions may apply; see official policy [here](https://www.kegg.jp/kegg/rest/)).

**Step 1**: Classify sequencing data to **KO** entries, bridged by a UniRef-to-KO mapping file:

```bash
woltka classify \
  --input  input_dir \
  --coords coords.txt.xz \
  --map    uniref/uniref.map.xz \
  --map    kegg/ko.map.xz \
  --rank   ko \
  --output ko.tsv
```

**Step 2**: Use the script [kegg_query.py](https://github.com/qiyunzhu/utils/blob/main/kegg_query.py) to build higher hierarchies of the KO's in the profile:

```bash
python kegg_query.py ko.tsv
```

Be patient, as the KEGG server limits the number of queries per time to 10 (see [policy](https://www.kegg.jp/kegg/rest/keggapi.html#list)).

This will generate multiple mapping files in the current directory. The filenames are self-explanatory. For examples: `ko-to-reaction.txt` is a mapping of KOs to [reactions](https://www.genome.jp/kegg/reaction/) (R), `ko-to-module.txt` is a mapping to [modules](https://www.genome.jp/kegg/module.html) (M), `ko-to-pathway.txt` is to [pathways](https://www.genome.jp/kegg/pathway.html), etc.

**Step 3**: Use Woltka's [**collapse**](collapse.md) command to convert KO's to higher-level classification units. This command supports many-to-many mapping, because that's the nature of the relationships between functional units (e.g., genes vs pathways). For example, the following command will generate a profile of **reactions**:

```bash
woltka tools collapse -i ko.tsv -m kegg/ko-to-reaction.txt -o reaction.tsv
```

Once a reaction profile is generated, several mapping files help to find more information about the reactions, such as `reaction_name.txt`, `reaction_equation.txt` and `reaction_definition.txt`. Meanwhile several other one-to-many mapping files, such as `reaction-to-compound.txt`, `reaction-to-module.txt`, `reaction-to-rclass.txt` help to explore the data toward other classification levels. All are enabled by the `collapse` command.

**Step 4**: You may want to explore the completeness of individual reactions, modules and pathways. This can be achieved using Woltka's [**coverage**](coverage.md) command. For example:

```bash
woltka tools coverage -i reaction.tsv -m kegg/module-to-reaction.txt -o module.cov.tsv
```

Differently from the last command, the **coverage** command generates a table, where cell values indicate the percentage of reactions required by each module found in each sample.

**Step 5**: You may also consider **stratifying** the functional profile by microbiome components. [See details](stratify.md).
