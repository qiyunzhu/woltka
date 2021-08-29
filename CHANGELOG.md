# Change Log

## Version 0.1.3 (08/27/2021)

### Changed
- Migrated from Travis CI to GitHub Actions (#127).
- Made `--map-as-rank` default when only mapping file(s) are provided (#132).
- Renamed `--normalize|-z` as `--frac|-f` (#128).
- Modified core algorithm which slightly improved performance (#124).

### Added
- Added `tool normalize` command, with multiple features (#124).
- Added the feature to collapse a stratified table (#126).
- Created an WoL FTP server, and added link to it (#118).
- Added an WoL standard operating procedure (`wolsop.sh`) and documentation (#116).
- Added first [citation](https://www.biorxiv.org/content/10.1101/2021.04.04.438427v1.abstract) of Woltka (#111).
- Added protocols for Bowtie2 / SHOGUN and Fastp (#121).
- Added discussion about mapping uniqueness (#131).

### Fixed
- Fixed free-rank classification subject not found issue (#120).
- Corrected paths to example files and directories (#117).


## Version 0.1.2 (03/31/2021)

### Changed
- Updated Qiita documentation (#107).
- Renamed "gOTU" with "OGU" (#104).

### Added
- Published at PyPI. Can be installed by `pip install woltka` (#108).
- Added instructions for using MetaCyc and KEGG (#99, #101).
- Added `tool collapse` command, which supports one-to-many classification (#99).

### Fixed
- Fixed Handling of zero length alignment (#105).


## Version 0.1.1 (02/17/2021)

### Added
- First official release.

