"""
Script to merge subindices search results into unified table

Output looks like this:

    S1  S2  S3  S4
q1  0.7     0.7 
q2      0.7     0.7
q3          0.7

Creates a sparse representation of binary
prsence/absence status. Missing values indicate
that query did not reach k-mer ratio during search.
"""

import os
import sys
import pandas as pd

def merge_sparse_dataframes(dataframes: list):
    base_df = pd.read_csv(dataframes[0], sep="\t", index_col=0).astype('Sparse[float64, nan]')
    for this_frame in range(1, len(dataframes)):
        tmp_df = pd.read_csv(dataframes[this_frame], sep="\t", index_col=0).astype('Sparse[float64, nan]')
        base_df  = pd.concat([base_df, tmp_df], axis=1)
    return base_df

merged_df = merge_sparse_dataframes(snakemake.input.indices)
merged_df.to_csv(snakemake.output.combined_indices, sep='\t', index=True)
#merged_df.to_hdf(snakemake.output.combined_indices_bin, index=True, key='merged_df')
