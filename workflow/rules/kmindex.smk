
rule ntcard:
    """
Estimate k-mer cardinality of each sample to be stored in index.
"""
    input:
        fastq=get_ntcard_fastq,
    output:
        histo="index/ntcard/{sample}.hist",
    log:
        "index/logs/ntcard/{sample}_ntcard.log",
    conda:
        "../envs/kmindex.yaml"
    container:
        "docker://quay.io/biocontainers/ntcard:1.2.2--pl5321hdcf5f25_4"
    threads: 1
    params:
        kmer_size=config["indexing"]["kmer_size"],
        prefix=lambda wildcards, output: output.histo.rstrip(".hist"),
    shell:
        "ntcard "
        "--kmer={params.kmer_size} "
        "--pref={params.prefix} "
        "{input.fastq} && "
        "mv index/ntcard/{wildcards.sample}_k*.hist index/ntcard/{wildcards.sample}.hist"


rule parse_ntcard:
    """
Parse ntcard k-mer histogram to get F0 and f1 counts.
"""
    input:
        histo="index/ntcard/{sample}.hist",
    output:
        parsed_histo="index/ntcard/{sample}.ntcard",
    log:
        "index/logs/ntcard/{sample}_ntcard_parsing.log",
    conda:
        "../envs/python3.yaml"
    container:
        "docker:://python:3.10.17-alpine3.22"
    threads: 1
    message:
        "Collecting k-mer cardinalities of all samples in indexing cohort"
    script:
        "../scripts/parse_ntcard.py"


rule gather_ntcard:
    """
Collect all k-mer cardinality counts, combine into one list and sort by largest number of unique k-mers.

"""
    input:
        histo=expand("index/ntcard/{sample}.ntcard", sample=samples.bin_id),
    output:
        kmer_all="index/ntcard/experiments.ntcard.txt",
        kmer_all_sorted="index/ntcard/experiments.ntcard.sorted.txt",
    log:
        "index/logs/ntcard/gather_ntcard.log",
    container:
        "docker://debian:bookworm-slim"
    threads: 1
    message:
        "Sorting k-mer cardinalities in descending order"
    shell:
        """
        printf 'experiment\\tF0\\tf1\\tnum_kmers\\n' > {output.kmer_all}
        cat {input.histo} >> {output.kmer_all}
        sort -nr -k 4 {output.kmer_all} > {output.kmer_all_sorted}
        """


rule estimate_bf_size:
    """
Estimate theoretical optimal bloom filter size.
"""
    input:
        kmer_all_experiments=rules.gather_ntcard.output.kmer_all_sorted,
    output:
        bloom_filter_size="index/kmindex/bloom_filter_size.txt",
    log:
        "index/logs/kmindex/bf_size_estimation.log",
    conda:
        "../envs/python2.yaml"
    container:
        "docker://python:2.7.18-alpine3.11"
    threads: 1
    params:
        bf_size_exe=workflow.source_path("../scripts/simple_bf_size_estimate.py"),
        fpr=float(config["indexing"]["fpr"]) * 100,
    message:
        "Estimating optimal Bloom filter size"
    shell:
        """
        largest_experiment="$(grep -v '#' {input.kmer_all_experiments} | head -n 1 | cut -f 4)"
        bf_size="$(python2 {params.bf_size_exe} $largest_experiment {params.fpr}% | grep -v '#' | head -n 1 | cut -f 4)"
        echo -n "$bf_size" > {output.bloom_filter_size}
        """


rule gather_fastq_kmtricks:
    input:
        fastq=get_ntcard_fastq,
    output:
        temp("index/kmindex/tmp/{sample}.txt"),
    log:
        "index/kmindex/{sample}_gather_fastq_input.log",
    container:
        "docker://debian:bookworm-slim"
    params:
        formatted_input=lambda wildcards, input: f"{wildcards.sample} : {' ; '.join(input.fastq)}",
    shell:
        'printf "{params.formatted_input}\n" > {output}'


rule write_kmtricks_fof:
    input:
        expand(
            "index/kmindex/tmp/{sample}.txt", sample=samples.bin_id.unique().tolist()
        ),
    output:
        fof="index/kmindex/samples.txt",
    log:
        "index/kmindex/create_fof.log",
    container:
        "docker://debian:bookworm-slim"
    shell:
        "cat {input} > {output.fof}"


rule kmtricks:
    """
Run kmtricks pipeline to extract and store k-mers as BFs.
"""
    input:
        bf_size=rules.estimate_bf_size.output.bloom_filter_size,
        sample_sheet=rules.write_kmtricks_fof.output.fof,
    output:
        kmtricks_index=directory("index/kmindex/kmtricks"),
    log:
        "index/kmindex/kmtricks_build.log",
    conda:
        "../envs/kmindex.yaml"
    container:
        "docker://tlemane/kmindex:0.5.2"
    threads: 16
    resources:
        mem_mb=20000,
    params:
        kmer_size=int(config["indexing"]["kmer_size"]),
        index_mode=(
            "hash:bfc:bin"
            if config["indexing"]["quantitative_index"]
            else "hash:bf:bin"
        ),
        abundance_classes=(
            f"-nb-cell 1000000 --bitw {config['indexing'].get('abundance_classes',2)}"
            if config["indexing"]["quantitative_index"]
            else ""
        ),
    shell:
        "test -d {output.kmtricks_index} && rmdir {output.kmtricks_index} ;"
        "bf_size=$(cat {input.bf_size}) && "
        "kmtricks "
        "pipeline "
        "--file {input.sample_sheet} "
        "--run-dir {output.kmtricks_index} "
        "--kmer-size {params.kmer_size} "
        "--hard-min 1 "
        "--mode {params.index_mode} "
        "--soft-min 2 "
        "--share-min 1 "
        "--bloom-size ${{bf_size}} "
        "--threads {threads} "
        "--minimizer-size 10 "
        "--nb-partitions 0 "
        "--cpr "
        "{params.abundance_classes} 2>&1 | tee {log}"


rule kmindex:
    """
Register kmtricks BF matrix into kmindex directory.
"""
    input:
        kmtricks_index=rules.kmtricks.output.kmtricks_index,
    output:
        kmindex_index=directory("index/kmindex/global_index"),
    log:
        "index/kmindex/kmindex_register.log",
    conda:
        "../envs/kmindex.yaml"
    container:
        "docker://tlemane/kmindex:0.5.2"
    threads: 1
    resources:
        mem_mb=10000,
    params:
        index_name=config["indexing"].get("index_name", "samples"),
    shell:
        """
        kmindex register -i {output.kmindex_index} -n {params.index_name} -p {input.kmtricks_index} 2>&1 | tee {log}
        """
