import os.path as path


rule query_raptor:
    """
Query raptor index with sequences in fasta file.
"""
    input:
        index=get_index,
        query_fasta=config["query"]["query_fasta"],
    output:
        search_results="query/raptor/{subindex}/search.tsv",
    log:
        "query/logs/raptor/{subindex}_search.log",
    conda:
        "../envs/raptor.yaml"
    container:
        "docker://quay.io/biocontainers/raptor:3.0.1--h6dccd9a_2"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input: get_memory_raptor(wildcards, input),
    params:
        theta=config["query"]["kmer_ratio"],
    shell:
        "raptor "
        "search "
        "--threshold {params.theta} "
        "--index {input.index} "
        "--query {input.query_fasta} "
        "--output {output.search_results} &> {log}"


rule kmindex_query:
    """
Query kmindex index with sequences in fasta file.
"""
    input:
        index=get_index,
        query_fasta=config["query"]["query_fasta"],
    output:
        search_results="query/kmindex/{subindex}/search.tsv",
    log:
        "query/logs/kmindex/{subindex}_search.log",
    conda:
        "../envs/kmindex.yaml"
    container:
        "docker://tlemane/kmindex:0.5.2"
    threads: 16
    params:
        output_dir=lambda wildcards, output: path.join(
            path.dirname(output.search_results), "search"
        ),
        findere_z=(
            f"--zvalue {config['query']['findere_z']}"
            if config["query"]["findere_z"] in [1, 2, 3, 4, 5, 6]
            else ""
        ),
    shell:
        "kmindex "
        "query "
        "-i {input.index} "
        "--fastx {input.query_fasta} "
        "--threads {threads} "
        "{params.findere_z} "
        "--aggregate "
        "--fast "
        "--format matrix "
        "--output {params.output_dir} &> {log}; "
        "mv {params.output_dir}/samples.tsv {output.search_results}"


rule parse_raptor_subindex_search:
    """
Parse Raptor search results into tabular format
"""
    input:
        search_results=rules.query_raptor.output.search_results,
    output:
        parsed_search_results="query/raptor/{subindex}/parsed_search.tsv.gz",
    log:
        "query/logs/raptor/{subindex}_parse_search.log",
    conda:
        "../envs/python3.yaml"
    threads: 1
    params:
        kmer_ratio=config["query"]["kmer_ratio"],
        exe=workflow.source_path("../scripts/parse_kmer_search.py"),
        index_mapping=lambda wildcards, input: index_struct[wildcards.subindex].get(
            "sample_mapping", ""
        ),
        uncompressed_file=lambda wildcards, output: output.parsed_search_results.rstrip(
            ".gz"
        ),
    shell:
        "python {params.exe} "
        "--search-results {input.search_results} "
        "--output {params.uncompressed_file} "
        "--method raptor "
        "--raptor-sample-mapping {params.index_mapping} "
        "--kmer-ratio {params.kmer_ratio} &> {log} && gzip -n {params.uncompressed_file} "


rule parse_kmindex_subindex_search:
    """
Parse Kmindex search results into tabular format
"""
    input:
        search_results=rules.kmindex_query.output.search_results,
    output:
        parsed_search_results="query/kmindex/{subindex}/parsed_search.tsv.gz",
    log:
        "query/logs/kmindex/{subindex}_parse_search.log",
    conda:
        "../envs/python3.yaml"
    threads: 1
    params:
        kmer_ratio=config["query"]["kmer_ratio"],
        exe=workflow.source_path("../scripts/parse_kmer_search.py"),
        uncompressed_file=lambda wildcards, output: output.parsed_search_results.rstrip(
            ".gz"
        ),
    shell:
        "python {params.exe} "
        "--search-results {input.search_results} "
        "--output {params.uncompressed_file} "
        "--method kmindex "
        "--kmer-ratio {params.kmer_ratio} &> {log} && gzip -n {params.uncompressed_file} "


rule gather_raptor_subindex_results:
    """
Combine raptor search results of sub-indices into unified table
"""
    input:
        indices=get_subindex_results_raptor,
    output:
        combined_indices_bin="query/raptor/search.parquet",
    log:
        "query/logs/raptor/merge_search.log",
    conda:
        "../envs/python3.yaml"
    threads: 1
    resources:
        mem_mb=20000,
    params:
        run_dir=lambda wildcards, output: os.path.dirname(output.combined_indices_bin),
    script:
        "../scripts/combined_kmer_search.py"


rule gather_kmindex_subindex_results:
    """
Combine kmindex search results of sub-indices into unified table
"""
    input:
        indices=get_subindex_results_kmindex,
    output:
        combined_indices_bin="query/kmindex/search.parquet",
    log:
        "query/logs/kmindex/merge_search.log",
    conda:
        "../envs/python3.yaml"
    threads: 1
    resources:
        mem_mb=20000,
    params:
        run_dir=lambda wildcards, output: os.path.dirname(output.combined_indices_bin),
    script:
        "../scripts/combined_kmer_search.py"
