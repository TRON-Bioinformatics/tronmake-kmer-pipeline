rule cuttlefish:
    input:
        fastq = get_ntcard_fastq()
    output:
        dbg = 'index/cuttlefish/{sample}_cdbg.fa'
        json = 'index/cuttlefish/{sample}_cdbg.json'
    threads: 8
    resources:
        mem_mb = 3000
    params:
        prefix = lambda wildcards, output:
            os.path.splitext(output.dbg)[0],
        input_csv = lambda wildcards, input:
            ','.join(input.fastq),
        kmer_size = int(config['indexing']['kmer_size']),
        cutoff = int(config['indexing']['cutoff'])
    conda:
        '../envs/cuttlefish.yaml'
    log:
        'index/logs/{sample}_cuttlefish.log'
    shell:
        'cuttlefish '
        '--read '
        '--kmer-len {params.kmer_size} '
        '--cutoff {params.cutoff} '
        '--seq={params.input_csv} '
        '--threads {threads} '
        '--output {params.prefix}'
        

