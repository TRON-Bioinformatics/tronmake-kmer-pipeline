rule bam2fastq:
    """
When BAM input is provided, separate reads into fastq files.
"""
    input:
        bam=get_bam_input,
    output:
        reads="index/prepare_input/{sample}/reads.fastq.gz",
    log:
        "index/logs/bam2fastq/{sample}_bam2fastq.log",
    conda:
        "../envs/samtools.yaml"
    container:
        "docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
    threads: 6
    shell:
        """
        exec 2> {log}
        # Check if multiple bam files are provided
        bam_files=({{{input.bam}}})
        no_bam_files=${{#bam_files[@]}}
        if [ $no_bam_files -gt 1 ]
        then
            samtools cat --threads 2 {input.bam} | samtools collate --threads 2 -u -O - | samtools fastq --threads 2 -o {output.reads} -0 /dev/null -n
        else
            samtools collate --threads 2 -u -O {input.bam} | samtools fastq --threads 2 -o {output.reads} -0 /dev/null -n
        fi
        """
