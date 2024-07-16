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
