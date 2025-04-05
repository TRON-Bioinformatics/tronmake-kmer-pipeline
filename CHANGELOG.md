# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2025-04-05

### Added

### Changed

* Updated pipeline to snakemake v9

* Added resource defaults for rules in workflow profile

* Updates to CI pipeline

* Updates to conda env yaml files

### Fixed

- Improved danymic memory allocation of raptor rules

### Known issues

### Removed

## [2.0.0] - 2024-08-08

### Added

* Added minimal WDL workflow for indexing in Terra cloud 

* Add expression classes for kmindex abundance indexing 

* Added support for BAM input e.g. TCGA and GTEx samples

* Added meta index support to query different indices in parallel

* Combined search results of subindices into detection matrix 

### Changed

* Ported pipeline to snakemake v8 

* Changed output of query format to compressed binary file (Parquet)

### Fixed

- Adjusted threads in raptor rules 

### Known issues

### Removed

