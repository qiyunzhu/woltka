# Change Log

## Version 0.1.6-dev

### Changed
- Changed default output subject coverage (`--outcov`) coordinates into BED-like (0-based, exclusive end). The output can be directly parsed by programs like bedtools. Also added support for GFF-like and other custom formats, as controled by paramter `--outcov-fmt` ([#204](https://github.com/qiyunzhu/woltka/pull/204) and [#205](https://github.com/qiyunzhu/woltka/pull/205)).


## Version 0.1.6 (2/22/2024)

### Changed
- Improved performance moderately ([#192](https://github.com/qiyunzhu/woltka/pull/192)).
- Parameter `--chunk` is now the number of unique query sequences instead of the number of lines ([#192](https://github.com/qiyunzhu/woltka/pull/192)).
- Updated GitHub Actions workflow.

### Added
- Added parameter `-x|--exclude`, which will exclude query sequences that are mapped to given reference sequences (such as host genome, spike-in, vector, etc.) ([#192](https://github.com/qiyunzhu/woltka/pull/192)).
- Added support for interleaved paired-end SAM files ([#191](https://github.com/qiyunzhu/woltka/pull/191)).
- Added native support for PAF file format ([#182](https://github.com/qiyunzhu/woltka/pull/182)).

### Fixed
- Updated hyperlinks in documentation.


## Version 0.1.5 (12/25/2022)

### Changed
- Modified the command-line interface. Sub-commands under the `tools` menu were raised to the main menu. For example, `woltka tools collapse` now becomes `woltka collapse` ([#176](https://github.com/qiyunzhu/woltka/pull/176)). The `tools` menus is still kept for backward compatibility ([#177](https://github.com/qiyunzhu/woltka/pull/177)).
- Improved efficiency of handling BIOM tables ([#171](https://github.com/qiyunzhu/woltka/pull/171), [#175](https://github.com/qiyunzhu/woltka/pull/175)). This significantly reduced the memory consumption and runtime of the `collapse` command.
- The `--trim-sub` parameter now takes underscore (`_`) as the default separator, if not explicitly specified ([#173](https://github.com/qiyunzhu/woltka/pull/173)).
- Updated installation protocol. Now Woltka can be Conda-installed without explicit pre-installation of biom-format ([#174](https://github.com/qiyunzhu/woltka/pull/174)).

### Added
- Upgraded the `collapse` command. Now it can collapse using internal hierarchies of feature IDs, such as genome-gene pairs, taxonomic lineages, and EC numbers. This upgrade enables stratified taxonomic/functional analysis of coord-matching ORF tables, without explicitly stratifying them during the classification step (which otherwise is time and space-consuming). The output should be identical to the that of the old method ([#173](https://github.com/qiyunzhu/woltka/pull/173)).
- Added a feature to the `normalize` command, which can extract gene lengths from a gene coordinates file. This enables convenient length-based normalization of an existing gene (ORF) profile. For example, it can convert raw frequencies into RPK ([#179](https://github.com/qiyunzhu/woltka/pull/179)).
- Extended description of stratification and collapsing protocols in the documentation.

### Fixed
- Fixed a bug in parsing sample ID lists ([#164](https://github.com/qiyunzhu/woltka/pull/164)).


## Version 0.1.4 (04/27/2022)

### Changed
- Used a single integer to store gene and read information in the ordinal mapper; use bitwise operations to parse the information. This significantly reduced memory consumption ([#142](https://github.com/qiyunzhu/woltka/pull/142)).
- Improved alignment file processing efficiency ([#153](https://github.com/qiyunzhu/woltka/pull/153)).
- Updated installation protocol. Currently a single `pip install` command should suffice ([#145](https://github.com/qiyunzhu/woltka/pull/145)).
- Updated GitHub Actions test from Python 3.6 to Python 3.8. The program itself should continue to support Python 3.6+ ([#158](https://github.com/qiyunzhu/woltka/pull/158)).

### Added
- Added a naive algorithm for read-gene matching when the number of reads is small. This improves speed ([#148](https://github.com/qiyunzhu/woltka/pull/148)).
- Added support for stdin as input. This lets the program take a variety of previously unsupported alignment formats ([#155](https://github.com/qiyunzhu/woltka/pull/155)).
- Added or updated several pieces of documentation, such as a RefSeq tutorial ([#157](https://github.com/qiyunzhu/woltka/pull/157)).
- Experimentally added a much accelerated ordinal mapper, powered by NumPy and Numba (branch `numba`) ([#152](https://github.com/qiyunzhu/woltka/pull/152)).


## Version 0.1.3 (08/27/2021)

### Changed
- Migrated from Travis CI to GitHub Actions ([#127](https://github.com/qiyunzhu/woltka/pull/127)).
- Made `--map-as-rank` default when only mapping file(s) are provided ([#132](https://github.com/qiyunzhu/woltka/pull/132)).
- Renamed `--normalize|-z` as `--frac|-f` ([#128](https://github.com/qiyunzhu/woltka/pull/128)).
- Modified core algorithm which slightly improved performance ([#124](https://github.com/qiyunzhu/woltka/pull/124)).

### Added
- Added `tool normalize` command, with multiple features ([#124](https://github.com/qiyunzhu/woltka/pull/124)).
- Added the feature to collapse a stratified table ([#126](https://github.com/qiyunzhu/woltka/pull/126)).
- Created an WoL FTP server, and added link to it ([#118](https://github.com/qiyunzhu/woltka/pull/118)).
- Added an WoL standard operating procedure (`wolsop.sh`) and documentation ([#116](https://github.com/qiyunzhu/woltka/pull/116)).
- Added first [citation](https://www.biorxiv.org/content/10.1101/2021.04.04.438427v1.abstract) of Woltka ([#111](https://github.com/qiyunzhu/woltka/pull/111)).
- Added protocols for Bowtie2 / SHOGUN and Fastp ([#121](https://github.com/qiyunzhu/woltka/pull/121)).
- Added discussion about mapping uniqueness ([#131](https://github.com/qiyunzhu/woltka/pull/131)).

### Fixed
- Fixed free-rank classification subject not found issue ([#120](https://github.com/qiyunzhu/woltka/pull/120)).
- Corrected paths to example files and directories ([#117](https://github.com/qiyunzhu/woltka/pull/117)).


## Version 0.1.2 (03/31/2021)

### Changed
- Updated Qiita documentation ([#107](https://github.com/qiyunzhu/woltka/pull/107)).
- Renamed "gOTU" with "OGU" ([#104](https://github.com/qiyunzhu/woltka/pull/104)).

### Added
- Published at PyPI. Can be installed by `pip install woltka` ([#108](https://github.com/qiyunzhu/woltka/pull/108)).
- Added instructions for using MetaCyc and KEGG ([#99](https://github.com/qiyunzhu/woltka/pull/99), [#101](https://github.com/qiyunzhu/woltka/pull/101)).
- Added `tool collapse` command, which supports one-to-many classification ([#99](https://github.com/qiyunzhu/woltka/pull/99)).

### Fixed
- Fixed Handling of zero length alignment ([#105](https://github.com/qiyunzhu/woltka/pull/105)).


## Version 0.1.1 (02/17/2021)

### Added
- First official release.
