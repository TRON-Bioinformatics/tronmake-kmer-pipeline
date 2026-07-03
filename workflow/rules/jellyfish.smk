checkpoint split_fasta:
    input:
        query_fasta=config["query"]["query_fasta"],
    output:
        directory("query/jellyfish/split_fasta"),
    log:
        "query/logs/jellyfish/split_fasta.log",
    container:
        "docker://busybox:1.36.1-musl"
    threads: 1
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
        query_fasta="query/jellyfish/split_fasta/{cts}.fasta",
        index=get_index,
    output:
        search_results=temp("query/jellyfish/{subindex}/{cts}_query.tsv"),
    log:
        "query/logs/jellyfish/{subindex}_{cts}_query.log",
    conda:
        "../envs/jellyfish.yaml"
    container:
        "docker://quay.io/biocontainers/kmer-jellyfish"
    threads: 1
    resources:
        mem_mb=8000,
    shell:
        "jellyfish query -s {input.query_fasta} -o {output.search_results} {input.index} 2> {log}"


rule jellyfish_parse:
    input:
        query="query/jellyfish/{subindex}/{cts}_query.tsv",
    output:
        parsed_result=temp("query/jellyfish/{subindex}/{cts}_parsed.tsv"),
    log:
        "query/logs/jellyfish/{subindex}_{cts}_query_parse.log",
    container:
        "docker://busybox:1.36.1-musl"
    threads: 1
    shell:
        """
        awk -v OFS='\\t' '
            {{
                print "{wildcards.cts}", $1, $2
            }}' {input.query} > {output.parsed_result} 2> {log}
        """


rule combine_jellyfish:
    input:
        aggregate_jellyfish_input,
    output:
        "query/jellyfish/{subindex}/quantitative_search.tsv",
    log:
        "query/logs/jellyfish/{subindex}_combine.log",
    container:
        "docker://busybox:1.36.1-musl"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"


rule jellyfish_index:
    input:
        fastq=get_ntcard_fastq,
    output:
        jf="index/jellyfish/{sample}.jf",
    log:
        "index/logs/jellyfish/{sample}_jf_build.log",
    conda:
        "../envs/jellyfish.yaml"
    container:
        "docker://quay.io/biocontainers/kmer-jellyfish"
    threads: 2
    resources:
        mem_mb=8000,
    params:
        kmer_size=int(config["indexing"]["kmer_size"]),
        workdir=lambda wildcards, output: os.path.dirname(output.jf),
    shell:
        """
        # Create FIFO as replacement for /dev/fd0
        ! test -f {params.workdir}/jf_fifo && mkfifo {params.workdir}/jf_fifo

        # Start gunzip in background and write to FIFO 
        zcat {input.fastq} > {params.workdir}/jf_fifo &
        # Save PID of gunzip command
        pid=$!

        jellyfish count \\
            {params.workdir}/jf_fifo \\
            -m {params.kmer_size} \\
            -s 1G \\
            -t {threads} \\
            -o {output.jf} \\
            &> {log}

        # Wait until zcat is finished before deleting the FIFO
        wait $pid
        rm -f {params.workdir}/jf_fifo
        """
