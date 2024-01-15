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
        