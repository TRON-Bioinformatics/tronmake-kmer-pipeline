import pandas as pd
import yaml

def validate_config(config):
    """
    Perform some validity checks on config parameters of pipeline.
    This ensures that selected options fit with workflow. There might be better options.
    """
    if config["modus"]["query"] and config["modus"]["indexing"]:
        raise ValueError("Can not run pipeline in indexing and query modus")
    if not config["modus"]["query"] and not config["modus"]["indexing"]:
        raise ValueError("Neither query not index mode activated. One run mode required")
    

    # Check user selected k-mer method is supported
    assert config['query']['method'] in ['raptor', 'kmindex'], "Selected method not supported"

    if config["modus"]["query"]:
        assert config['query']['index'] is not None, \
            "k-mer index required for search"
        assert config['query']['kmer_ratio'] is not None, \
            "k-mer ratio required for search"

    if config["modus"]["indexing"]:
        assert config["indexing"]['samples'] is not None, \
            "Sample table required for indexing"
        assert int(config["indexing"]['kmer_size']) in list(range(19,32)), \
            "k-mer size not supported"

def check_index_struct(index_struct):
    """
    Dummy method. Requires proper checking of index struct
    """
    return True

def verbose_logs(verbose: bool = True):
    pass

def get_final_output(samples: pd.DataFrame, index_struct: dict):
    """
    Populate final output for target rule. Final output is determined
    based on selected run mode.
    """
    final_output = []
    if config['modus']['query']:
        for index_id, index_properties in index_struct.items():
            method = index_properties.get('method', None)
            if method is None:
                raise ValueError(f"k-mer method not specified for index: {index_id}")
            final_output.append(f'query/{method}/{index_id}/search.txt')
    
    elif config['modus']['indexing']:
        method = config['indexing']['method']
        match method:
            case "raptor":
                final_output.append(
                    'index/raptor/raptor.index'
                )
                final_output.append(
                    'index/raptor/index_mapping.txt'
                )
            case "kmindex":
                final_output.append(
                    'index/kmindex/global_index'
                )
            case "cuttlefish":
                final_output.extend(
                    expand('index/cuttlefish/{sample}/{sample}_cdbg.fa.gz', sample = samples.bin_id)
                )
            case _:
                sys.exit(1)
    else:
        print("Neither indexing or querying selected.")
        sys.exit(1)

    return final_output

def read_sample_sheet(file):
    """
    Read sample sheet with fastq files for indexing
    """
    file_content = []
    with open(file, "r") as file_handle:
        for line in file_handle:
            elements = line.rstrip().split('\t')
            bin_id = elements[0].rstrip()
            fastq = elements[1].rstrip()
            file_content.append({'bin_id': bin_id, 'fastq': fastq})
    samples = pd.DataFrame(file_content)
    return samples

def read_index_struct(file):
    """
    Read meta index of different raptor indices
    """
    with open(file, 'r') as file_handle:
        index_struct = yaml.safe_load(file_handle)
    # ToDo add sanity checks if meta index file is correctly formatted
    if not check_index_struct(index_struct):
        sys.exit(1)
    return index_struct

def get_ntcard_fastq(wildcards):
    sample = samples.query('bin_id == @wildcards.sample')
    fastq = sample.get('fastq')
    try: 
        fastq = fastq.item()
    except AttributeError as error:
        return ''
    # Our sample sheet format allows the following delimiters ","
    fastq = fastq.split(',')
    return fastq

def get_number_of_bin(samples):
    """
    Get number of samples in input sample sheet
    """
    return samples.shape[0]
    

def get_memory_raptor(wildcards, input):
    """
    Calculate memory required for raptor queries.
    Memory consumption is approximately the index size on disk. We add some additional 
    memory to prevent the job from failing
    """

    memory=max(input.size_mb * 1.1, input.size_mb)
    return memory

def get_memory_raptor_build(wildcards):
    """
    Calculate memory required for raptor build step.
    Memory consumption is growing approximately linearly with number of bins to be indexed.
    The scaling factor was determined by indexing multiple batches of sequencing data to assess
    scalability. We add some additional memory to prevent the job from failing.
    scaling_factor: ~0.15 
    """
    samples_to_index = get_number_of_bin(samples)
    estimated_memory = round(samples_to_index * 0.15 * 1024)
    estimated_memory = max(estimated_memory, 150000)
    return estimated_memory

def get_index(wildcards):
    index_to_query = index_struct.get(wildcards.subindex)
    path = index_to_query.get('path', '')
    return path
