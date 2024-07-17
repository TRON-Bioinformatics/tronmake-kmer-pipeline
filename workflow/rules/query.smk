import os.path as path


rule query_raptor:
    """
    Query raptor index with sequences in fasta file.
    """
    input:
        index = get_index,
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/raptor/{subindex}/search.txt'
    params:
        theta = config['query']['kmer_ratio']
    threads: 1
    resources:
        mem_mb = lambda wildcards, input: get_memory_raptor(wildcards, input)
    conda:
        '../envs/raptor.yaml'
    log: 'query/raptor/{subindex}_search.log'
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
        search_results = 'query/kmindex/{subindex}/search.txt'
    params:
        output_dir = lambda wildcards, output: path.join(path.dirname(output.search_results), "search"),
        findere_z = f"--zvalue {config['query']['findere_z']}" \
            if config['query']['findere_z'] in [1,2,3,4,5,6] else ''
    threads: 16
    conda:
        '../envs/kmindex.yaml'
    log: 'query/kmindex/{subindex}_search.log'
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
    input:
        search_results = rules.query_raptor.output.search_results
    output:
        parsed_search_results = query/raptor/{subindex}/parsed_search.tsv
    threads: 1
    params:
        kmer_ratio = config['query']['kmer_ratio'],
        exe = workflow.source_path("../scripts/parse_kmer_search.py"),
        index_mapping = lambda wildcards, input: index_struct[wildcards.subindex].get('index_mapping', '')
    shell:
        'python {params.exe} '
        '--search-results {input.search_results} '
        '--output {output.parsed_search_results} '
        '--method raptor '
        '--raptor-sample-mapping {params.index_mapping}'

rule parse_kmindex_subindex_search:
    input:
        search_results = rules.kmindex_query.output.search_results
    output:
        parsed_search_results = 'query/kmindex/{subindex}/parsed_search.tsv',
        kmer_ratio = config['query']['kmer_ratio'],
        exe = workflow.source_path("../scripts/parse_kmer_search.py"),
    threads: 1
    shell:
        'python {params.exe} '
        '--search-results {input.search_results} '
        '--output {output.parsed_search_results} '
        '--method kmindex '
        '--kmindex-cutoff {prams.kmer_ratio}'

rule gather_raptor_subindex_results:
    input:
        indices = get_subindices_raptor
    output:
        combined_indices = query/raptor/search.tsv
    shell:
        '''
        cat {input.indices} > query/raptor/concat.tmp
        sort -k1,1 < query/raptor/concat.tmp > {output.combined_indices}
        rm query/raptor/concat.tmp
        '''


