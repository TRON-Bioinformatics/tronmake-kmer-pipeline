# TronMake k-mer pipeline

<!-- badges: start -->

[![Release](https://img.shields.io/badge/release-v0.0.1-blue?style=flat)](https://gitlab.rlp.net/tron/kmer_pipeline)
[![Snakemake](https://img.shields.io/badge/snakemake-7.31.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

<!-- badges: end -->


The TronMake k-mer pipeline is part of the software solution k4neo for indexing large collections of healthy and tumor
RNA-seq datasets to demonstrate tumor-specificity of (neo)antigen candidates.

The workflow implements the best practice workflows of the state-of-the-art k-mer indexing methods [Raptor](https://github.com/seqan/raptor)
and [kmindex](https://github.com/tlemane/kmindex).

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




## How to run it

In general we recommend running the workflow from the k4neo python package. However, you can also
execute it manually to create k-mer indices and to execute queries against it. 

Download the project and run as follows. Modify the [config file](config/config.yaml) to fit your needs
and execute the workflow (1). You can also change the parameters on the command line (2).

```bash

# (1) Run pipeline using custom config file
snakemake --use-conda -j 24

# Create k-mer index with raptor
snakemake --config modus={indexing:true} indexing={samples:/path/to/your/samples, method:raptor}

```

### Input tables

The table with fastq files expects two tab-separated columns with a header. Multiple FASTQs can be provided separated by commas.

| bin_id   | fastq                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz                              |
| sample_2 | /path/to/sample_2.fastq.gz,/path/to/sample_2_2.fastq.gz |

## References

* Mehringer, S., Seiler, E., Droop, F. et al. Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries. Genome Biol 24, 131 (2023). https://doi.org/10.1186/s13059-023-02971-4

* Téo Lemane, Nolan Lezzoche, Julien Lecubin, Eric Pelletier, Magali Lescot, Rayan Chikhi, Pierre Peterlongo kmindex and ORA: indexing and real-time user-friendly queries in terabytes-sized complex genomic datasets. bioRxiv 2023.05.31.543043; doi: https://doi.org/10.1101/2023.05.31.543043 