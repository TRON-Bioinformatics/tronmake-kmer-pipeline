rule cuttlefish:
    """
    Extract compacted DBG from sequencing samples
    """
    input:
        fastq = get_ntcard_fastq
    output:
        dbg = temp('index/cuttlefish/{sample}/{sample}_cdbg.fa'),
        json = 'index/cuttlefish/{sample}/{sample}_cdbg.json'
    threads: 2
    resources:
        mem_mb = 4000
    params:
        prefix = lambda wildcards, output:
            os.path.splitext(output.dbg)[0],
        work_dir = lambda wildcards, output: os.path.dirname(output.dbg),
        input_csv = lambda wildcards, input:
            ','.join(input.fastq),
        kmer_size = int(config['indexing']['kmer_size']),
        cutoff = int(config['indexing']['cutoff'])
    conda:
        '../envs/cuttlefish.yaml'
    container:
        'docker://quay.io/biocontainers/cuttlefish:2.2.0--h6a68c12_2'
    log:
        'index/logs/{sample}_cuttlefish.log'
    shell:
        'cuttlefish '
        'build '
        '--read '
        '--kmer-len {params.kmer_size} '
        '--work-dir {params.work_dir} '
        '--cutoff {params.cutoff} '
        '--seq={params.input_csv} '
        '--threads {threads} '
        '--output {params.prefix} &> {log}'

rule compress_cDBG:
    input:
        dbg = rules.cuttlefish.output.dbg
    output:
        compress_dbg = 'index/cuttlefish/{sample}/{sample}_cdbg.fa.gz'
    threads: 1
    container: 'docker://busybox:1.36.1-musl'
    shell:
        '''
        gzip {input.dbg}
        '''
        

