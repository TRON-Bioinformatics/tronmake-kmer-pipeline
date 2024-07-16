# TronMake k-mer pipeline

<!-- badges: start -->

[![Release](https://img.shields.io/badge/release-v1.0.0-blue?style=flat)](https://gitlab.rlp.net/tron/kmer_pipeline)
[![Snakemake](https://img.shields.io/badge/snakemake-8.16.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

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



## How to run it

In general we recommend running the workflow from the k4neo python package. However, you can also
execute it manually to create k-mer indices and to execute queries against it. 

Download the project and run as follows. Modify the [config file](config/config.yaml) to fit your needs
and execute the workflow. Moreover you need to provide an index manifest file in **query** mode describing the
indices to query. Have a look at the [example manifest](tests/configs/raptor_index.yaml)

We also provide a python wrapper to execute the workflow and configure the pipeline based on the arguments passed.

### Wrapper execution (recommended)

```bash

# Create k-mer index with raptor
python tronmake-kmer.py \
    index \
    --samples examples/samples.tsv \
    --kmer 21 \
    --method raptor \
    --fpr 0.05 \
    --workdir /path/to/index \
    --slurm

# Create k-mer index with kmindex
python tronmake-kmer.py \
    index \
    --samples examples/samples.tsv \
    --kmer 21 \
    --method kmindex \
    --fpr 0.05 \
    --workdir /path/to/index \
    --slurm

# Query k-mer indices in manifest file 
python tronmake-kmer.py \
    query \
    --index-manifest tests/configs/raptor_index.yaml \
    --fasta examples/query.fasta \
    --detection-ratio 0.7 \
    --workdir /path/to/query \
    --slurm


```



### Traditional execution

```bash

# (1) Run pipeline using custom config file
snakemake --use-conda -j 24

# Create k-mer index with raptor
snakemake --config modus={indexing:true} indexing={samples:/path/to/your/samples, method:raptor}

```

### Input

#### Index mode

The table with fastq files expects two tab-separated columns with a header. Multiple FASTQs can be provided separated by commas.

| bin_id   | fastq                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz                              |
| sample_2 | /path/to/sample_2.fastq.gz,/path/to/sample_2_2.fastq.gz |


#### Query mode

The index manifest file for indexing expects a YAML file describing the indices, their location and optionally a sample-index mapping.

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


## References

* Mehringer, S., Seiler, E., Droop, F. et al. Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries. Genome Biol 24, 131 (2023). https://doi.org/10.1186/s13059-023-02971-4

* Téo Lemane, Nolan Lezzoche, Julien Lecubin, Eric Pelletier, Magali Lescot, Rayan Chikhi, Pierre Peterlongo kmindex and ORA: indexing and real-time user-friendly queries in terabytes-sized complex genomic datasets. bioRxiv 2023.05.31.543043; doi: https://doi.org/10.1101/2023.05.31.543043  

* Jamshed Khan, Rob Patro, Cuttlefish: fast, parallel and low-memory compaction of de Bruijn graphs from large-scale genome collections, Bioinformatics, Volume 37, Issue Supplement_1, July 2021, Pages i177–i186, https://doi.org/10.1093/bioinformatics/btab309  