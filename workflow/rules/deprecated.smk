rule reindeer_query:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/reindeer/reindeer_search.txt'
    params:
        output_dir = lambda wildcards, output: path.join(path.dirname(output.search_results), "search"),
        reindeer_exe = config['query']['reindeer_exe'],
        theta = config['query']['kmer_ratio']
    threads: 1
    log: 'query/reindeer/search.log'
    shell:
        '{params.reindeer_exe} '
        '--query '
        '-P  {params.theta} '
        '-l {input.index} '
        '-q {input.query_fasta} '
        '-o {params.output_dir} &> {log}; '
        'mv'

rule cobs_query:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/cobs/cobs_search.txt'
    params:
        theta = config['query']['kmer_ratio']
    threads: 1
    conda:
        '../envs/cobs.yaml'
    log: 'query/cobs/search.log'
    shell:
        'cobs '
        'query '
        '--index {input.index} '
        '--file {input.query_fasta} '
        '--threshold {params.theta} '
        '> {output.search_results} 2> {log}'