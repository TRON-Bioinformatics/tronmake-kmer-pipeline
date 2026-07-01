# tronmake k-mer pipeline

<!-- badges: start -->

[![Snakemake](https://img.shields.io/badge/snakemake-9.23.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
![Pytest](https://img.shields.io/badge/pytest-%23ffffff.svg?style=flat-square&logo=pytest&logoColor=2f9fe3)
[![CI](https://github.com/TRON-Bioinformatics/tronmake-kmer-pipeline/actions/workflows/ci.yaml/badge.svg)](https://github.com/TRON-Bioinformatics/tronmake-kmer-pipeline/actions/workflows/ci.yaml/badge.svg)
[![Release](https://img.shields.io/badge/release-v3.0.0-blue?style=flat)](https://github.com/TRON-Bioinformatics/tronmake-kmer-pipeline/releases/tag/v3.0.0)
[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](https://opensource.org/licenses/MIT)

<!-- badges: end -->


The tronmake k-mer pipeline is a Snakemake workflow designed to index large collections of sequencing data using different k-mer datastructures.

The pipeline implements best-practice workflows for several state-of-the-art k-mer indexing methods: [Raptor](https://github.com/seqan/raptor), [kmindex](https://github.com/tlemane/kmindex), compacted De Bruijn Graph (cDBG) construction via [cuttlefish](https://github.com/COMBINE-lab/cuttlefish) and counting bloom filter construction via [Jellyfish](https://github.com/gmarcais/Jellyfish)

Depending on the selected method, the workflow follows different processing steps:

**Raptor**
* Extract minimisers from raw FASTA files (using a window of $k+4$).
* Generate a Hierarchical Interleaved Bloom Filter (HIBF) layout.
* Construct the HIBF index based on the specified $w, k$ schema.

**Kmindex**
* Estimate k-mer cardinality for individual samples.
* Determine the optimal Bloom filter size based on estimated cardinality.
* Extract k-mers using Kmtricks to create either presence/absence or quantitative indices.
* Construct a global meta-index using Kmindex.

**Cuttlefish**
* Extract k-mers and construct compacted De Bruijn Graphs (cDBG).
* Compress maximal unitigs for efficient storage and querying.

**Jellyfish**
* Construct k-mer counting filters from individual samples.
* Search these filters for query k-mers.
* Parse and annotate each k-mer count using the input FASTA files.

## Usage

To run it, download the project, create the conda environment and adapt the [config file](tests/configs/kmindex.yaml) to fit your requirements. 

```
conda env create -p conda_env -f environment.yaml
```

When running in **query** mode, you must provide an index manifest file describing the indices to be queried. See the [example manifest](tests\index_manifest\raptor_index_manifest.yaml) for details.

### Input

#### Index Mode

The input is a tab-separated sample sheet **without** a header. Multiple files per sample can be provided, separated by commas.

**Option 1: FASTQ only**

If only FASTQ files are used, the table should contain two columns: `bin_id` and `fastq` file path.

| bin_id   | fastq                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz                              |
| sample_2 | /path/to/sample_2.fastq.gz,/path/to/sample_2_2.fastq.gz |

**Option 2: Mixed FASTQ and BAM**

If the input includes reads in (u)BAM format, the table must contain three columns: `bin_id`, `fastq` (or `bam`) file path, and `file_type`.

| bin_id   | fastq/bam                                               | file_type |
|:--------:|:-------------------------------------------------------:|:---------:|
| sample_1 | /path/to/sample_1.fastq.gz                              | fastq     |
| sample_2 | /path/to/sample_2.bam,/path/to/sample_2_2.bam           | bam       |
| sample_3 | /path/to/sample_3.fastq.gz,/path/to/sample_3_2.fastq.gz | fastq     |

#### Query Mode

For querying, an index manifest file (YAML) is required to describe the indices, their locations, and optionally a sample-index mapping for parsing Raptor results. The query sequences themselves are specified in the general [config file](config/config.yaml).

Example manifest:
```yaml
raptor_test_index:
  samples: 1
  path: 'examples/raptor.index'
  sample_mapping: 'examples/index_mapping.txt'
  method: raptor

kmindex_test_index:
  samples: 1
  path: 'index/kmindex/global_index'
  method: kmindex

jellyfish_test_index:
  path: 'examples/jellyfish/IT_N_103.jf'
  method: jellyfish
```

Using this manifest, the pipeline will search for query sequences across all listed indices.


### Execution

#### Snakemake Command Line

```bash

snakemake \
    --directory <output dir> \
    [--software-deployment-method conda \]
    [--software-deployment-method apptainer \]
    --configfile </path/to/myconfig> \
    [--conda-prefix </path/to/conda/env/location>]

```

* `directory`: Specifies where the query/index results are stored.
* `configfile`: The path to the config file.
* `software-deployment-method`: Currently conda and apptainer are supported and tested.
* `--conda-prefix` (optional): Where should the conda environments be stored

### Output

#### Index mode

The output is contained in the `index` directory inside the folder specified with the option `--directory`.
The resulting indices files are separated by method.

For kmindex the index directory would look like the following structure:
```
index/
├── kmindex
│   ├── bloom_filter_size.txt
│   ├── global_index
│   ├── kmindex_register.log
│   ├── kmtricks
│   ├── kmtricks_build.log
│   └── samples.txt
├── ntcard
│   ├── 1.hist
│   ├── 1.ntcard
│   ├── 2.hist
│   ├── 2.ntcard
│   ├── experiments.ntcard.sorted.txt
│   └── experiments.ntcard.txt
└── prepare_input
    ├── 1
    └── 2
```

* `index/kmindex/bloom_filter_size.txt` : Estimated optiomal Bloom filter size based on k-mer cardinality of largest sample.
* `index/kmindex/global_index` : Global kmindex index
* `index/kmindex/kmtricks` : Kmtricks count/binary Bloom filters (required by global index)
* `index/kmindex/samples.txt` : File of files sample sheet required by kmtricks
* `index/kmindex/ntcard` : K-mer cardinality of indexed samples.
* `index/kmindex/prepare_input` : Contains fastq files of BAM samples.

For Raptor the index directory would look like the following structure:
```
index/
├── prepare_input
│   ├── 1
│   └── 2
└── raptor
    ├── build.log
    ├── hibf_binning.layout
    ├── index_mapping.txt
    ├── layout.log
    ├── minimiser
    └── raptor.index
```

* `index/raptor/raptor.index` : Raptor HIBF index
* `index/raptor/minimiser` : Raptor winnowing minimiser files
* `index/raptor/index.mapping` : Sample mapping. This file maps raptor bin paths to sample names in the original sample sheet (required for parsing results of query)
* `index/prepare_input` : Contains fastq files of BAM samples.


#### Query mode

The output is contained in the `query` directory inside the folder specified with the option `--directory`.
The resulting search results are organized by sub-indices defined in the k-mer manifest file. 
A sample directory for a given index search contains the subdirectories of subindices and aggregated search results:

```
query/
└── kmindex
    ├── kmindex_test_index
    │   ├── parsed_search.tsv.gz
    │   ├── parse_search.log
    │   ├── search
    │   ├── search.log
    │   └── search.tsv
    ├── merge_search.log
    └── search.parquet
```


Aggregated search results are provided as compressed binary file. You can use R or python to read the detection matrix file.

* `query/<method>/search.parquet`: Contains aggregated search results over all subindices as Apache parquet files.

Quantitative annotation of search results are provided for each JellyFish index as simple TSV file.

* `query/jellyfish/<index_name>/quantitative_search.tsv`

## Cloud Execution (WDL)

For cloud-based execution, a WDL workflow is provided in the `WDL/` directory. This allows for scalable indexing of large datasets using platforms like Terra or Cromwell, handling tasks such as BAM to FASTQ conversion and distributed minimiser extraction before building the final Raptor HIBF index. Other methods are currently not supported.


## References

* Mehringer, S., Seiler, E., Droop, F. et al. Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries. Genome Biol 24, 131 (2023). https://doi.org/10.1186/s13059-023-02971-4

* Téo Lemane, Nolan Lezzoche, Julien Lecubin, Eric Pelletier, Magali Lescot, Rayan Chikhi, Pierre Peterlongo kmindex and ORA: indexing and real-time user-friendly queries in terabytes-sized complex genomic datasets. bioRxiv 2023.05.31.543043; doi: https://doi.org/10.1101/2023.05.31.543043  

* Jamshed Khan, Rob Patro, Cuttlefish: fast, parallel and low-memory compaction of de Bruijn graphs from large-scale genome collections, Bioinformatics, Volume 37, Issue Supplement_1, July 2021, Pages i177–i186, https://doi.org/10.1093/bioinformatics/btab309  

* Marçais, G., & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics, 27(6), 764-770.

* Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., van Dyken, P. C., Tomkins-Tinch, C. H., ... & Köster, J. (2025). Sustainable data analysis with Snakemake. F1000Research, 10, 33.
