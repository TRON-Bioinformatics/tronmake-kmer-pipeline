import os.path as path


rule query_raptor:
    """
    Query raptor index with sequences in fasta file.
    """
    input:
        index = get_index,
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/raptor/{subindex}/search.tsv'
    params:
        theta = config['query']['kmer_ratio']
    threads: 1
    resources:
        mem_mb = lambda wildcards, input: get_memory_raptor(wildcards, input)
    conda:
        '../envs/raptor.yaml'
    log: 'query/raptor/{subindex}/search.log'
    shell:
        'raptor '
        'search '
        '--threshold {params.theta} '
        '--index {input.index} '
        '--query {input.query_fasta} '
        '--output {output.search_results} &>{log}'

rule kmindex_query:
    """
    Query kmindex index with sequences in fasta file.
    """
    input:
        index = get_index,
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/kmindex/{subindex}/search.tsv'
    params:
        output_dir = lambda wildcards, output: path.join(path.dirname(output.search_results), "search"),
        findere_z = f"--zvalue {config['query']['findere_z']}" \
            if config['query']['findere_z'] in [1,2,3,4,5,6] else ''
    threads: 16
    conda:
        '../envs/kmindex.yaml'
    log: 'query/kmindex/{subindex}/search.log'
    shell:
        'kmindex '
        'query '
        '-i {input.index} '
        '--fastx {input.query_fasta} '
        '--threads {threads} '
        '{params.findere_z} '
        '--aggregate '
        '--fast '
        '--format matrix '
        '--output {params.output_dir} &> {log}; '
        'mv {params.output_dir}/samples.tsv {output.search_results}'

rule parse_raptor_subindex_search:
    """
    Parse Raptor search results into tabular format
    """
    input:
        search_results = rules.query_raptor.output.search_results
    output:
        parsed_search_results = 'query/raptor/{subindex}/parsed_search.tsv.gz'
    threads: 1
    params:
        kmer_ratio = config['query']['kmer_ratio'],
        exe = workflow.source_path("../scripts/parse_kmer_search.py"),
        index_mapping = lambda wildcards, input: index_struct[wildcards.subindex].get('sample_mapping', ''),
        uncompressed_file = lambda wildcards, output: output.parsed_search_results.rstrip('.gz')
    log: 'query/raptor/{subindex}/parse_search.log'
    conda:
        '../envs/python3.yaml'
    shell:
        'python {params.exe} '
        '--search-results {input.search_results} '
        '--output {params.uncompressed_file} '
        '--method raptor '
        '--raptor-sample-mapping {params.index_mapping} '
        '--kmer-ratio {params.kmer_ratio} && gzip -n {params.uncompressed_file} '

rule parse_kmindex_subindex_search:
    """
    Parse Kmindex search results into tabular format
    """
    input:
        search_results = rules.kmindex_query.output.search_results
    output:
        parsed_search_results = 'query/kmindex/{subindex}/parsed_search.tsv.gz'
    threads: 1
    params:    
        kmer_ratio = config['query']['kmer_ratio'],
        exe = workflow.source_path("../scripts/parse_kmer_search.py"),
        uncompressed_file = lambda wildcards, output: output.parsed_search_results.rstrip('.gz')
    threads: 1
    conda:
        '../envs/python3.yaml'
    log: 'query/kmindex/{subindex}/parse_search.log'
    shell:
        'python {params.exe} '
        '--search-results {input.search_results} '
        '--output {params.uncompressed_file} '
        '--method kmindex '
        '--kmer-ratio {params.kmer_ratio} && gzip -n {params.uncompressed_file} '

rule gather_raptor_subindex_results:
    """
    Combine raptor search results of sub-indices into unified table
    """
    input:
        indices = get_subindex_results_raptor
    output:
        combined_indices_bin = 'query/raptor/search.parquet'
    params:
        run_dir = lambda wildcards, output: os.path.dirname(output.combined_indices_bin)
    conda:
        '../envs/python3.yaml'
    log: 'query/raptor/merge_search.log'
    threads: 1
    resources:
        mem_mb = 20000
    script:
        '../scripts/combined_kmer_search.py'

rule gather_kmindex_subindex_results:
    """
    Combine kmindex search results of sub-indices into unified table
    """
    input:
        indices = get_subindex_results_kmindex
    output:
        combined_indices_bin = 'query/kmindex/search.parquet'
    params:
        run_dir = lambda wildcards, output: os.path.dirname(output.combined_indices_bin)
    threads: 1
    resources:
        mem_mb = 20000
    conda:
        '../envs/python3.yaml'
    log: 'query/kmindex/merge_search.log'
    script:
        '../scripts/combined_kmer_search.py'
