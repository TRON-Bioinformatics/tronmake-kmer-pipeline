
rule query_raptor:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/raptor/raptor_search.txt'
    params:
        theta = config['query']['kmer_ratio']
    threads: 1
    conda:
        '../envs/raptor.yaml'
    log: 'query/raptor/search.log'
    shell:
        'raptor '
        'search '
        '--threshold {params.theta} '
        '--index {input.index} '
        '--query {input.query_fasta} '
        '--output {output.search_results} &>{log}'

rule kmindex_query:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/kmindex/kmindex_search.txt'
    params:
        output_dir = lambda wildcards, output: path.join(path.dirname(output.search_results), "search")
    threads: 1
    conda:
        '../envs/kmindex.yaml'
    log: 'query/kmindex/search.log'
    shell:
        'kmindex '
        'query '
        '-i {input.index} '
        '--fastx {input.query_fasta} '
        '--fast '
        '--format matrix '
        '--output {params.output_dir} &> {log}; '
        'mv {params.output_dir}/samples.tsv {output.search_results}'

