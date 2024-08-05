# TronMake k-mer pipeline

<!-- badges: start -->

[![Release](https://gitlab.rlp.net/tron/tronmake-kmer-pipeline/-/badges/release.svg)](https://gitlab.rlp.net/tron/tronmake-kmer-pipeline/-/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-8.16.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![pipeline status](https://gitlab.rlp.net/tron/tronmake-kmer-pipeline/badges/develop/pipeline.svg)](https://gitlab.rlp.net/tron/tronmake-kmer-pipeline/commits/master)


<!-- badges: end -->


The TronMake k-mer pipeline is part of the software solution k4neo for indexing large collections of healthy and tumor
RNA-seq datasets to demonstrate tumor-specificity of (neo)antigen candidates.

The workflow implements the best practice workflows of the state-of-the-art k-mer indexing methods [Raptor](https://github.com/seqan/raptor)
and [kmindex](https://github.com/tlemane/kmindex) as well as compacted De Bruijn Graphs (cDBG) construction with [cuttlefish](https://github.com/COMBINE-lab/cuttlefish).

Depending on the selected k-mer method the workflow consist of different steps.

**Raptor**

* Extracting minimiser from raw fasta files (`k+4`)
* Creating a HIBF layout 
* Creating HIBF index with `w,k` schema

**Kmindex**

* Estimate k-mer cardinality of individual samples
* Estimate optimal Bloom filter size based on k-mer cardinality
* Extracting k-mers with Kmtricks (presence/absence or quantitative index)
* Creation of global meta-index with Kmindex

**Cuttlefish**

* Extract k-mers and construct compacted De Bruijn Graphs (cDBG)
* Compress maximal unitigs



## Usage

In general we recommend running the workflow from the k4neo python package. However, you can also
execute it manually to create k-mer indices and to execute queries against it. 

Download the project and run as follows. Modify the [config file](config/config.yaml) to fit your needs
and execute the workflow. Moreover you need to provide an index manifest file in **query** mode describing the
indices to query. Have a look at the [example manifest](tests/configs/raptor_index.yaml)

We also provide a python wrapper (python >= 3.10 supported) to execute the workflow and configure the pipeline based on the arguments passed.

### Input

#### Index mode

The table with fastq files expects two (three) tab-separated columns without a header. Multiple FASTQs can be provided separated by commas.

| bin_id   | fastq                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz                              |
| sample_2 | /path/to/sample_2.fastq.gz,/path/to/sample_2_2.fastq.gz |

If input also contains reads in (u)BAM format, the table with input files expects three tab-separated columns without a header. Multiple FASTQs and BAMs can be provided separated by commas.

| bin_id   | fastq                                                   | file_type |
|:--------:|:-------------------------------------------------------:|:---------:|
| sample_1 | /path/to/sample_1.fastq.gz                              | fastq     |
| sample_2 | /path/to/sample_2.bam,/path/to/sample_2_2.bam           | bam       |
| sample_3 | /path/to/sample_3.fastq.gz,/path/to/sample_3_2.fastq.gz | bam       |

#### Query mode

The index manifest file for querying expects a YAML file describing the indices, their location and optionally a sample-index mapping required for parsing Raptor results.
The query sequences are specified in the general [config file](config/config.yaml).

```
raptor_test_index:
  samples: 1
  path: 'examples/raptor.index'
  sample_mapping: 'examples/index_mapping.txt'
  method: raptor

kmindex_test_index:
  samples: 1
  path: 'index/kmindex/global_index'
  method: kmindex
```

With this manifest the query sequences would be searched in the k-mer indices of Raptor and kmindex.


### Execution

#### SnakeMake command line

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
* `software-deployment-method`: Currently conda and singularity is supported and tested.
* `--conda-prefix` (optional): Where should the conda environments be stored

#### Wrapper execution (recommended)

```bash

# Create a k-mer index with raptor
python tronmake-kmer.py \
    index \
    --samples examples/samples.tsv \
    --kmer 21 \
    --method raptor \
    --fpr 0.05 \
    --workdir /path/to/index \
    --slurm

# Create a k-mer index with kmindex
python tronmake-kmer.py \
    index \
    --samples examples/samples.tsv \
    --kmer 21 \
    --method kmindex \
    --fpr 0.05 \
    --workdir /path/to/index \
    --slurm

# Query k-mer indices defined in manifest file for sequences in fasta file
python tronmake-kmer.py \
    query \
    --index-manifest tests/index_manifest/raptor_index_manifest.yaml \
    --fasta examples/query.fasta \
    --detection-ratio 0.7 \
    --workdir /path/to/query \
    --slurm

```

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
* `index/kmindex/prepare_input` : Contains fastq files of BAM samples.


#### Query mode

The output is contained in the `query` directory inside the folder specified with the option `--directory`.
The resulting search files are separated by sub-indices defined in the k-mer manifest file. 
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

* `query/<method>/search.parquet`: Containes aggregated search results over all subindices as Apache parquet files.


## References

* Mehringer, S., Seiler, E., Droop, F. et al. Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries. Genome Biol 24, 131 (2023). https://doi.org/10.1186/s13059-023-02971-4

* Téo Lemane, Nolan Lezzoche, Julien Lecubin, Eric Pelletier, Magali Lescot, Rayan Chikhi, Pierre Peterlongo kmindex and ORA: indexing and real-time user-friendly queries in terabytes-sized complex genomic datasets. bioRxiv 2023.05.31.543043; doi: https://doi.org/10.1101/2023.05.31.543043  

* Jamshed Khan, Rob Patro, Cuttlefish: fast, parallel and low-memory compaction of de Bruijn graphs from large-scale genome collections, Bioinformatics, Volume 37, Issue Supplement_1, July 2021, Pages i177–i186, https://doi.org/10.1093/bioinformatics/btab309  