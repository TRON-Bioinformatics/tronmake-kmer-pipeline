import pandas as pd


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

def verbose_logs(verbose: bool = True):
    pass

def get_final_output():
    """
    Populate final output for target rule. Final output is determined
    based on selected run mode.
    """
    final_output = []

    if config['modus']['query']:
        method = config['query']['method']
        final_output.append(
            f'query/{method}/{method}_search.txt'
        )
    
    elif config['modus']['indexing']:
        method = config['query']['method']
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
            case _:
                sys.exit(1)
    else:
        print("Neither indexing or querying selected.")
        sys.exit(1)

    return final_output

def read_sample_sheet(file):
    """
    Read sample sheet as input
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

def get_ntcard_fastq(wildcards):
    sample = samples.query('bin_id == @wildcards.sample')
    fastq = sample.get('fastq')
    try: 
        fastq = fastq.item()
    except AttributeError as error:
        return ''
    # Our sample sheet format allows the following delimiters ";"
    fastq = fastq.split(',')
    return fastq

def get_memory_raptor(wildcards, input):
    """
    Calculate memory required for raptor queries.
    Memory consumption is approximately the index size on disk. We add some additional 
    memory to prevent thew job from failing
    """

    memory=max(input.index.size_mb * 1.1, input.index.size_mb)
    return memory
