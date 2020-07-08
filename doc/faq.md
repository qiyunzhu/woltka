# Frequently asked questions

## Input files

### Can Woltka parse compressed files?

Yes. All input files for Woltka (alignments and databases) can be supplied as compressed in gzip, bzip2 or xz formates. Woltka will automatically recognize and parse them.

## Hierarchies

### What if some subjects are not found in the classification system?

It is okay. Tolerating them and reporting them as "unassigned" is by design and not a bug in Woltka.


## Computing

### How many CPU cores can Woltka use?

Woltka works the best with two CPU cores (threads): one for decompression and the other for classification. This is automatically determined.