import os.path


rule raptor_prepare_per_sample:
    """
    Extract k+4 minimisers for each input sample
    """
    input:
        bin_fastq=get_ntcard_fastq,
    output:
        bin_minimiser="index/raptor/minimiser/{sample}/minimiser.list",
        idx_map="index/raptor/minimiser/{sample}/index_map.txt",
    params:
        sample_sheet=lambda wildcards, output: os.path.join(
            os.path.dirname(output.bin_minimiser), "sample_sheet.txt"
        ),
        output_dir=lambda wildcards, output: os.path.dirname(output.bin_minimiser),
        cut_off=int(config["indexing"]["cutoff"]),
        kmer_size=int(config["indexing"]["kmer_size"]),
        window=int(config["indexing"]["kmer_size"]) + 4,
    threads: 1
    conda:
        "../envs/raptor.yaml"
    container:
        "docker://quay.io/biocontainers/raptor:3.0.1--h6dccd9a_2"
    log:
        "index/logs/raptor/minimiser/{sample}_minimiser.log",
    message:
        "Extracting {params.kmer_size},{params.window} minimisers from sample {wildcards.sample}"
    shell:
        """
        exec 2> {log}
        echo {input.bin_fastq} > {params.sample_sheet}
        raptor prepare --threads {threads} --kmer {params.kmer_size} \\
        --window {params.window} --kmer-count-cutoff {params.cut_off} \\
        --input {params.sample_sheet} \\
        --output {params.output_dir} &> {log}
        paste {params.output_dir}/minimiser.list <(echo {wildcards.sample}) > {params.output_dir}/index_map.txt
        """


rule gather_raptor_minimisers:
    """
    Collect minimiser files of all samples to be included in the index.
    """
    input:
        minimisers=expand(
            "index/raptor/minimiser/{sample}/minimiser.list",
            sample=samples.bin_id.unique().tolist(),
        ),
    output:
        minimiser_list="index/raptor/minimiser/minimiser.list",
    threads: 1
    container:
        "docker://busybox:1.36.1-musl"
    log:
       "index/logs/raptor/minimiser/minimiser_gathering.log",
    shell:
        """
        exec 2> {log}
        cat {input.minimisers} > {output.minimiser_list}
        """


rule raptor_sample_mapping:
    input:
        idx_map=expand(
            "index/raptor/minimiser/{sample}/index_map.txt",
            sample=samples.bin_id.unique().tolist(),
        ),
    output:
        index_mapping="index/raptor/index_mapping.txt",
    threads: 1
    container:
        "docker://busybox:1.36.1-musl"
    log:
        "index/logs/raptor/index_mapping.log",
    shell:
        """
        exec 2> {log}
        {{ printf "minimiser_id\tsample_name\n" ; cat {input.idx_map} ; }} > {output.index_mapping}
        """


rule raptor_layout:
    """
    Create HIBF layout file from minimiser files.WS
    """
    input:
        minimiser_list=rules.gather_raptor_minimisers.output.minimiser_list,
    output:
        layout_file="index/raptor/hibf_binning.layout",
    params:
        fpr=float(config["indexing"]["fpr"]),
    conda:
        "../envs/raptor.yaml"
    container:
        "docker://quay.io/biocontainers/raptor:3.0.1--h6dccd9a_2"
    log:
        "index/logs/raptor/layout.log",
    threads: 1
    resources:
        mem_mb=500,
    message:
        "Determining HIBF index layout"
    shell:
        "raptor "
        "layout "
        "--input-file {input.minimiser_list} "
        "--false-positive-rate {params.fpr} "
        "--threads {threads} "
        "--output-filename {output.layout_file} &> {log}"


rule raptor_build:
    """
    Build Raptor HIBF index from layout and minimiser files.
    """
    input:
        layout_file=rules.raptor_layout.output.layout_file,
    output:
        hibf_index="index/raptor/raptor.index",
    conda:
        "../envs/raptor.yaml"
    container:
        "docker://quay.io/biocontainers/raptor:3.0.1--h6dccd9a_2"
    log:
        "index/logs/raptor/build.log",
    threads: 16
    resources:
        mem_mb=get_memory_raptor_build,
    message:
        "Building Raptor index"
    shell:
        "raptor "
        "build "
        "--input {input.layout_file} "
        "--threads {threads} "
        "--compressed "
        "--output {output.hibf_index} &> {log}"
