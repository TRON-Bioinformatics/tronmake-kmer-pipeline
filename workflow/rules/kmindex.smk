import os.path as path

rule ntcard:
    input:
        fastq = get_ntcard_fastq
    output:
        histo = "index/ntcard/{sample}.hist"
    params:
        kmer_size = config['indexing']['kmer_size'],
        prefix = lambda wildcards, output: output.histo.rstrip(".hist")
    conda:
        '../envs/kmindex.yaml'
    shell:
        'ntcard '
        '--kmer={params.kmer_size} '
        '--pref={params.prefix} '
        '{input.fastq} && '
        'mv index/ntcard/{wildcards.sample}_k*.hist index/ntcard/{wildcards.sample}.hist'

rule parse_ntcard:
    input:
        histo = "index/ntcard/{sample}.hist"
    output:
        parsed_histo = "index/ntcard/{sample}.ntcard"
    run:
        with open(input.histo, 'r') as file_handle, open(output.parsed_histo, 'w') as write_handle:
            f0, f1 = 0, 0
            for line in file_handle:
                elements = line.rstrip().split("\t")
                if elements[0] == "F0":
                    f0 = int(elements[1])
                elif elements[0] == "1":
                    f1 = int(elements[1])
                else:
                    continue
            result = f"{wildcards.sample}\t{f0}\t{f1}\t{f0-f1}\n"
            write_handle.write(result)

rule gather_ntcard:
    input:
        histo = expand("index/ntcard/{sample}.ntcard", sample = samples.bin_id)
    output:
        kmer_all = 'index/ntcard/experiments.ntcard.txt',
        kmer_all_sorted = 'index/ntcard/experiments.ntcard.sorted.txt'
    shell:
        '''
        printf 'experiment\\tF0\\tf1\\tnum_kmers\\n' > {output.kmer_all}
        cat {input.histo} >> {output.kmer_all}
        sort -nr -k 4 {output.kmer_all} > {output.kmer_all_sorted}
        '''

rule estimate_bf_size:
    input:
        kmer_all_experiments = rules.gather_ntcard.output.kmer_all_sorted
    params:
        bf_size_exe = '../scripts/simple_bf_size_estimate.py',
        fpr = float(config['indexing']['fpr']) * 100
    output:
        bloom_filter_size = 'index/kmindex/bloom_filter_size.txt'
    shell:
        '''
        largest_experiment="$(grep -v '#' {input.kmer_all_experiments} | head -n 1 | cut -f 4)"
        bf_size="$(python2 {params.bf_size_exe} $largest_experiment {params.fpr}% | grep -v '#' | head -n 1 | cut -f 4)"
        echo -n "$bf_size" > {output.bloom_filter_size}
        '''

rule write_kmtricks_fof:
    output:
        fof = 'index/kmindex/samples.txt'
    run:
        with open(output.fof, 'w') as file_handle:
            for sample in samples.itertuples():
                fastq = sample.fastq.replace(",", ";")
                file_handle.write(f'{sample.bin_id}:{fastq}\n')

rule kmtricks:
    input:
        bf_size = rules.estimate_bf_size.output.bloom_filter_size,
        sample_sheet = rules.write_kmtricks_fof.output.fof
    output:
        kmtricks_index = directory('index/kmindex/kmtricks')
    threads: 16
    params:
        kmer_size = int(config['indexing']['kmer_size']),
    conda:
        '../envs/kmindex.yaml'
    shell:
        'test -d {output.kmtricks_index} && rmdir {output.kmtricks_index} ;'
        'bf_size=$(cat {input.bf_size}) && '
        'kmtricks '
        'pipeline '
        '--file {input.sample_sheet} '
        '--run-dir {output.kmtricks_index} '
        '--kmer-size {params.kmer_size} '
        '--hard-min 1 '
        '--mode hash:bf:bin '
        '--soft-min 2 '
        '--share-min 1 '
        '--bloom-size ${{bf_size}} '
        '--threads {threads} '
        '--minimizer-size 10 '
        '--nb-partitions 0 '
        '--cpr'

rule kmindex:
    input:
        kmtricks_index = rules.kmtricks.output.kmtricks_index
    output:
        kmindex_index = directory('index/kmindex/global_index')
    conda:
        '../envs/kmindex.yaml'
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        '''
        # Delete automatically generated output dir
        kmindex register -i {output.kmindex_index} -n samples -p {input.kmtricks_index}
        '''

