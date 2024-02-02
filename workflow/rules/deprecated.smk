import os.path as path

rule bcalm2:
    """
    Transform collection of fastq files into compacted De Bruijn Graphs
    """
    input:
        fastq = get_fastq_for_bin
    params:
        kmer_size = config['indexing']['kmer_size'],
        cutoff = config['indexing']['cutoff']
    shadow: "shallow"
    output:
        cdbg = "index/bcalm2/graphs/{sample}.unitigs.fa.gz"
    shell:
        'bcalm '
        '-nb-cores {threads} '
        '-in {input.fastq} '
        '-kmer-size {params.kmer_size} '
        '-abundance-min {params.cutoff} '
        '-max-memory {resources.mem_mb} '
        '-out {wildcards.sample} && '
        'gzip {wildcards.sample}.unitigs.fa'

rule collect_graphs:
    input:
        cdbg = expand("index/bcalm2/graphs/{sample}.unitigs.fa.gz",
            sample = samples.bin_id)
    output:
        dbg_fof = "index/bcalm2/graphs_lst.txt"
    run:
        with open(output.dbg_fof, "w") as write_handle:
            for unitig in input.cdbg:
                write_handle.write(unitig + '\n')
        


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

rule convert_bin_to_raptor_fof:
    """Parse sample sheet to FOF

    k4neo supports the kmindex sheet format as input. However, Raptor uses
    it's own file of files format. This rule converts the sample sheet for raptor

    """
    input:
        sample_sheet = config['indexing']['samples']
    output:
        fof = "index/raptor/fof.txt",
        index_mapping = "index/raptor/index_mapping.txt"
    run:
        with open(input.sample_sheet, 'r') as file_handle, open(output.fof, 'w') as write_handle, open(output.index_mapping, 'w') as mapping_handle :

            mapping_handle.write("sample_name\tminimiser_id\n")
            for line in file_handle:
                elements = line.rstrip().split(' : ')
                fastq = elements[1].replace(";", " ")
                
                # Write FOF
                write_handle.write(fastq + "\n")
                # Write index mapping: sample_name: minimider_id
                # Raptor uses the basename of the first fastq file as bin identifier
                minimiser_id = os.path.basename(fastq.split(" ")[0]).rstrip(".minimiser")
                mapping_handle.write(f'{elements[0]}\t{minimiser_id}\n')

rule raptor_prepare:
    """
    Prepare FASTQ files for raptor indexing
    """
    input:
        sample_sheet = rules.convert_bin_to_raptor_fof.output.fof,
    output:
        minimiser_list = 'index/raptor/minimiser/minimiser.list'
    params:
        cut_off = int(config['indexing']['cutoff']),
        kmer_size = int(config['indexing']['kmer_size']),
        window = int(config['indexing']['kmer_size']) + 4,
        output_dir = lambda wildcards, output: os.path.dirname(output.minimiser_list)
    conda:
        '../envs/raptor.yaml'
    threads: 16
    resources:
        mem_mb = 50000
    log:
        'index/raptor/minimiser/minimiser_creation.log'
    shell:
        'raptor '
        'prepare '
        '--threads {threads} '
        '--kmer {params.kmer_size} '
        '--window {params.window} '
        '--kmer-count-cutoff {params.cut_off} '
        '--input {input.sample_sheet} '
        '--output {params.output_dir} &> {log}'