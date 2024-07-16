version 1.0

workflow raptor {
    meta {
        title: "Raptor"
        summary: "Index large collections of RNA-seq data with Raptor using HIBF index"
        authors: "Johannes Hausmann"
        version: "13"
    }
    input {
        Array[String] input_ids
        Array[File] bam_input_files
        Int kmer = 21
        Int minimiser_window
        Int recurrence_min = 2
        Int cpu
        Int cpu_pre_processing = 4
        Float fpr = 0.05
        Int hash = 2
        Boolean use_ssd = false
        Int preemptible = 3
    }
    parameter_meta {
        input_ids: "Unique sample identifier used as bin in Raptor index"
        bam_input_files: "Reads in BAM file"
        kmer: "K-mer size used to construct index"
        minimiser_window: "Size of minimiser window"
        recurrence_min: "Minimum recurrence of k-mer in sequencing reads to be included in index"
        cpu: "Number of CPUs to construct k-mer index"
        cpu_pre_processing: "Number of CPUs for preprocessing steps"
        fpr: "Theoretical FPR of Bloom filters in index"
        hash: "Number of hash functions in Bloom filter"
        preemptible: "Number of preemptible VM attempts before falling back to regular VM"
    }
    Int window = if !defined(minimiser_window) then kmer + 4 else minimiser_window
    String minimiser_directory = "raptor_minimiser"
    # Scatter bamtofastq over samples --> probably split Raptor and samtools to speed up things but then we have to pay storage cost
    scatter(idx in range(length(input_ids))) {
        call bam2fastq {
            input:
                input_id = input_ids[idx],
                bam_file = bam_input_files[idx],
                cpu = cpu_pre_processing,
                use_ssd = use_ssd,
                preemptible = preemptible
                
        }
        call raptor_prepare_per_sample {
            input:
                input_id = input_ids[idx],
                fastq_file = bam2fastq.fastq,
                minimiser_directory = minimiser_directory,
                cpu = cpu_pre_processing,
                kmer = kmer,
                window = window,
                recurrence_min = recurrence_min,
                use_ssd = use_ssd,
                preemptible = preemptible
        }

    }
    Array[File] minimiser_files = raptor_prepare_per_sample.minimiser
    Array[File] header_files = raptor_prepare_per_sample.header
    Array[File] fastqs = bam2fastq.fastq
    # Raptor HIBF index is between 2-5% of original input size in GB. We give it a bit more disk space (15%) to hold also minimiser files
    Int disk_size_index_build = ceil(size(fastqs, "GiB") * 0.15) + ceil(size(minimiser_files, "GiB"))

    call raptor_build {
        input:
            minimisers = minimiser_files,
            headers = header_files,
            fpr = fpr,
            hash = hash,
            cpu = cpu,
            number_of_bins = length(input_ids),
            disk_size_gb = disk_size_index_build,
            use_ssd = use_ssd,
            preemptible = 0
    }
    output {
        File layout = raptor_build.layout
        File raptor_index = raptor_build.raptor_index
    }
}

task raptor_prepare_per_sample {
    input {
        String input_id
        File fastq_file
        String minimiser_directory
        Int cpu
        Int kmer
        Int window
        Int recurrence_min
        Boolean use_ssd
        Int preemptible
        Int dynamic_disk_size = ceil(size(fastq_file,"GiB")) + 50
    }
    runtime {
        cpu: cpu
        memory: "20G"
        maxRetries: 1
        disks: "local-disk " + dynamic_disk_size + (if use_ssd then " SSD" else " HDD")
        docker: 'quay.io/biocontainers/raptor:3.0.1--h6dccd9a_2'
        preemptible: preemptible
    }
    output {
        File header = "~{minimiser_directory}/~{input_id}.header"
        File minimiser = "~{minimiser_directory}/~{input_id}.minimiser"
    }
    command <<<
        set -e
        echo ~{fastq_file} > ~{input_id}.lst
        raptor prepare \
            --input ~{input_id}.lst \
            --output ~{minimiser_directory} \
            --kmer ~{kmer} \
            --window ~{window} \
            --threads 4 \
            --kmer-count-cutoff ~{recurrence_min}
        rm ~{input_id}.lst
    >>>
}

task bam2fastq {
    input {
        String input_id
        File bam_file
        Int cpu
        Boolean use_ssd
        Int dynamic_disk_size = ceil(size(bam_file,"GiB"))*2 + 50
        Int preemptible
    }
    runtime {
        cpu: cpu
        memory: "10G"
        memory_gb: "10G"
        disks: "local-disk " + dynamic_disk_size + (if use_ssd then " SSD" else " HDD")
        docker: 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'
        preemptible: preemptible
        maxRetries: 3
    }
    output {
        File fastq = "~{input_id}.fastq.gz"
    }
    command <<<
        set -e
        samtools collate --threads 2 -u -O ~{bam_file} | \
            samtools fastq --threads 2 -o ~{input_id}.fastq.gz -0 /dev/null
    >>>
}

task raptor_build {
    input {
        Array[File] minimisers
        Array[File] headers
        Int hash
        Float fpr
        Int cpu
        Int number_of_bins
        Int disk_size_gb
        Boolean use_ssd
        Int memory_gb = ceil(number_of_bins * 0.15) + 50
        Int preemptible
    }
    output{
        File raptor_index = "raptor.index"
        File layout = "raptor.layout"
    }
    runtime {
        cpu: cpu
        memory: "~{memory_gb}G"
        memory_gb: "~{memory_gb}G"
        disks: "local-disk ~{disk_size_gb} " + (if use_ssd then "SSD" else "HDD")
        docker: 'quay.io/biocontainers/raptor:3.0.1--h6dccd9a_2'
        preemptible: preemptible
    }
    command <<<
        set -e
        printf "~{sep='\n' minimisers}" > minimiser.lst
        raptor layout \
            --input-file minimiser.lst \
            --false-positive-rate ~{fpr} \
            --num-hash-functions ~{hash} \
            --output-filename raptor.layout
        raptor build \
            --input raptor.layout \
            --output raptor.index \
            --fpr ~{fpr} \
            --threads ~{cpu} \
            --hash ~{hash} \
            --compressed
    >>>
}