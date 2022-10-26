# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.4]
### Fixed
- Fix transcript matrices not in output folder.
### Added
- output of merged bam optional.
- Repeat umap creation with different random states.
### Changed 
- Transcript counting Salmon on stringtie-generated transcriptome.
- Several perforance-related reforactorings including reductions in read write operations. 
- single_cell_sample_sheet is optional and kit options can be supplied as workflow parameters.


## [v0.1.3]
## Changed
- Better handling of sample_id conflicts in single_cell_sampkle_sheet and fastgingress.
- single_cell_sample_sheet is optional.
- Minor IO performance enhancements.
### Added
- kit options can be supplied from command line/config and applied to all samples. 

## [v0.1.2]
### Added
- Transcript x cell matrix output.

## [v0.1.1]
### Changed
- Combined gathering and splitting of fastqs into a single process.
- Use split2 for splitting fastqs.
- Remove unused kneeplot flags.
### Added
- check for identical sample_ids in single cell sample sheet and fastq data.

## [v0.1.0]
### Added
- First release. Port of [Sockeye](https://github.com/nanoporetech/sockeye) to Nextflow


