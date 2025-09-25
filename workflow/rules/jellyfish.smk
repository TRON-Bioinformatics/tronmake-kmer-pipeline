checkpoint split_fasta:
    input:
        query_fasta = config["query"]["query_fasta"],
    output:
        directory("query/jellyfish/split_fasta")
    container:
        "docker://busybox:1.36.1-musl"
    shell:
        """
        mkdir -p {output}
        awk '
            /^>/ {{
                f = substr($1, 2) ".fasta"
                file = "{output[0]}/" f
                print $0 > file
                next
            }}
            {{
                print $0 >> file
            }}' {input.query_fasta}
        """

rule jellyfish_query:
    input:
        query_fasta = "query/jellyfish/split_fasta/{cts}.fasta",
        index = get_index
    output:
        search_results = temp("query/jellyfish/{subindex}/{cts}.tsv")
    threads: 1
    resources:
        mem_mb = 8000
    conda:
        '../envs/jellyfish.yaml'
    container:
        'docker://quay.io/biocontainers/kmer-jellyfish'
    shell:
       "jellyfish query -s {input.query_fasta} -o {output.search_results} {input.index}"

rule jellyfish_parse:
    input:
        query = "query/jellyfish/{subindex}/{cts}.tsv"
    output:
        parsed_result = temp("query/jellyfish/{subindex}/{cts}_parsed.tsv")
    container:
        "docker://busybox:1.36.1-musl"
    shell:
        """
        awk -v OFS='\\t' '{{print "{wildcards.cts}", $1, $2}}' {input.query} > {output.parsed_result}
        """

def aggregate_input(wildcards):
    
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    
    return expand("query/jellyfish/{subindex}/{cts}_parsed.tsv",
           subindex=wildcards.subindex,
           cts=glob_wildcards(os.path.join(checkpoint_output, "{cts}.fasta")).cts)

rule combine_jellyfish:
    input:
        aggregate_input
    output:
        "query/jellyfish/{subindex}/quantitative_search.tsv"
    container:
        "docker://busybox:1.36.1-musl"
    shell:
        "cat {input} > {output}"


