checkpoint split_fasta:
    input:
        query_fasta = config["query"]["query_fasta"],
    output:
        directory("query/jellyfish/split_fasta")
    container:
        "docker://busybox:1.36.1-musl"
    threads: 1
    log:
        "query/logs/jellyfish/split_fasta.log",
    shell:
        """
        mkdir -p {output}
        awk '
            /^>/ {{
                f = substr($1, 2) ".fasta"
                file = "{output}/" f
                print $0 > file
                next
            }}
            {{
                print $0 >> file
            }}' {input.query_fasta} 2> {log}
        """

rule jellyfish_query:
    input:
        query_fasta = "query/jellyfish/split_fasta/{cts}.fasta",
        index = get_index
    output:
        search_results = temp("query/jellyfish/{subindex}/{cts}.tsv")
    threads: 1
    log:
        "query/logs/jellyfish/{subindex}_{cts}_query.log",
    resources:
        mem_mb = 8000
    conda:
        '../envs/jellyfish.yaml'
    container:
        'docker://quay.io/biocontainers/kmer-jellyfish'
    shell:
       "jellyfish query -s {input.query_fasta} -o {output.search_results} {input.index} 2> {log}"

rule jellyfish_parse:
    input:
        query = "query/jellyfish/{subindex}/{cts}.tsv"
    output:
        parsed_result = temp("query/jellyfish/{subindex}/{cts}_parsed.tsv")
    container:
        "docker://busybox:1.36.1-musl"
    threads: 1
    log:
        "query/logs/jellyfish/{subindex}_{cts}_query_parse.log",
    shell:
        """
        awk -v OFS='\\t' '
            {{
                print "{wildcards.cts}", $1, $2
            }}' {input.query} > {output.parsed_result} 2> {log}
        """

rule combine_jellyfish:
    input:
        aggregate_jellyfish_input
    output:
        "query/jellyfish/{subindex}/quantitative_search.tsv"
    container:
        "docker://busybox:1.36.1-musl"
    threads: 1
    log:
        "query/logs/jellyfish/{subindex}_combine.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule jellyfish_index:
    input:
        fastq = get_ntcard_fastq
    output:
        jf = "index/jellyfish/{sample}.jf"
    params:
        kmer_size = int(config["indexing"]["kmer_size"]),
    conda:
        '../envs/jellyfish.yaml'
    container:
        'docker://quay.io/biocontainers/kmer-jellyfish'
    threads: 2
    resources:
        mem_mb = 8000
    log:
        'index/logs/jellyfish/{sample}_jf_build.log'
    shell:
        'zcat {input.fastq} | '
        'jellyfish count '
        '/dev/stdin '
        '-m {params.kmer_size} '
        '-s 1G '
        '-t {threads} '
        '-o {output.jf} '
        '&> {log}'
    