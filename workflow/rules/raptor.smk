import os.path

rule raptor_prepare_per_sample:
    input:
        bin_fastq = get_ntcard_fastq
    output:
        bin_minimiser = 'index/raptor/minimiser/{sample}/minimiser.list'
    params:
        sample_sheet = 
            lambda wildcards, output: 
                os.path.join(os.path.dirname(output.bin_minimiser), 'sample_sheet.txt'),
        output_dir =
            lambda wildcards, output:
                os.path.dirname(output.bin_minimiser),
        cut_off = int(config['indexing']['cutoff']),
        kmer_size = int(config['indexing']['kmer_size']),
        window = int(config['indexing']['kmer_size']) + 4,
    threads: 4
    log: 'index/raptor/minimiser/{sample}/log'
    message: "Extracting {params.kmer_size},{params.window} minimisers from sample {wildcards.sample}"
    shell:
        '''
        echo {input.bin_fastq} > {params.sample_sheet}
        raptor prepare --threads {threads} --kmer {params.kmer_size} \\
        --window {params.window} --kmer-count-cutoff {params.cut_off} \\
        --input {params.sample_sheet} \\
        --output {params.output_dir} &> {log}
        '''

rule gather_raptor_minimisers:
    input:
        minimisers = expand('index/raptor/minimiser/{sample}/minimiser.list',
            sample=samples.bin_id.unique().tolist())
    output:
        minimiser_list= 'index/raptor/minimiser/minimiser.list'
    threads: 1
    shell:
        '''
        cat {input.minimisers} > {output.minimiser_list}
        '''

rule raptor_sample_mapping:
    """Parse sample sheet to FOF

    k4neo supports the kmindex sheet format as input. However, Raptor uses
    it's own file of files format. This rule converts the sample sheet for raptor

    """
    output:
        index_mapping = "index/raptor/index_mapping.txt"
    message: "Generating index-bin to sample mapping"
    run:
        with open(output.index_mapping, 'w') as mapping_handle :
            mapping_handle.write("sample_name\tminimiser_id\n")
            for line in samples.itertuples(index=False):
                # Write index mapping: sample_name: minimider_id
                # Raptor uses the basename of the first fastq file as bin identifier
                minimiser_id = os.path.basename(line.fastq.split(",")[0]).rstrip(".fastq.gz")
                mapping_handle.write(f'{line.bin_id}\t{minimiser_id}\n')

rule raptor_layout:
    """
    Create HIBF layout file from minimiser files
    """
    input:
        minimiser_list = rules.gather_raptor_minimisers.output.minimiser_list
    output:
        layout_file = "index/raptor/hibf_binning.layout"
    params:
        fpr = float(config['indexing']['fpr'])
    conda:
        '../envs/raptor.yaml'
    log: 'index/raptor/layout.log'
    threads: 16
    resources:
        mem_mb = 500
    message: "Determining HIBF index layout"
    shell:
        'raptor '
        'layout '
        '--input-file {input.minimiser_list} '
        '--false-positive-rate {params.fpr} '
        '--threads {threads} '
        '--output-filename {output.layout_file} &> {log}'

rule raptor_build:
    """
    Build Raptor HIBF index from layout and minimiser files
    """
    input:
        layout_file = rules.raptor_layout.output.layout_file
    output:
        hibf_index = "index/raptor/raptor.index"
    conda:
        '../envs/raptor.yaml'
    log: 'index/raptor/build.log'
    threads: 16
    resources:
        mem_mb=150000
    message: "Building Raptor index"
    shell:
        'raptor '
        'build '
        '--input {input.layout_file} '
        '--threads {threads} '
        '--compressed '
        '--output {output.hibf_index} &> {log}'
