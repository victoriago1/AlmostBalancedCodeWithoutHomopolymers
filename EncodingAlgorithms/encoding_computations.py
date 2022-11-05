from math import log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
sys.path.append('EncodingAlgorithms/')
import strand_requirements
import common
import RLL3andKnuth
import RLL3andKnuthwithSkips
import RLL3andKnuthbyBlocks

"""
A script for running the calculations of baseline and proposed solutions,
and saving them in separate and common CSV files.
"""

def compute_log2_count(strand_redundancy_df, method):
    """ Computes for every length of strands that hold both constraints,
        the log2 of the possible quaternary strands that can be representes by that length.
        (approximately equals to: the maximal original binary data length that it can represent) """
    start = 1
    log2_count = np.zeros(strand_requirements.MAX_n_quaternary + start + 1)

    print("Calculating log2 count for strands up to {} bits for {}".format(strand_requirements.MAX_n_binary, method))

    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary + start)):
        # calculate the maximal length that may represent a binary string of length n:
        last_original_length = strand_redundancy_df.loc[strand_redundancy_df['Final_length'] <= n,
                                                        'Original_length'].max()

        # when strand_redundancy_array[last_original_length, COL_FINAL_LENGTH] >= n,
        # the length n is sufficient to represent quaternary strands of length last_original_length,
        # meaning log2(4^last_original_length) = 2*last_original_length possible strands.
        log2_count[n] = 2*last_original_length

    log2_count = pd.DataFrame(log2_count)
    log2_count.index.name = "Final Quaternary Strand Length"
    log2_count.columns = [method + " log2 of possible strands"]
    log2_count.to_csv("EncodingAlgorithms/Results/" + method + "- Log2 Count.csv")
    return log2_count


def compute_quaternary_redundancy(temp_strand_redundancy_array, method):
    start = 1
    redundancy = np.zeros(strand_requirements.MAX_n_quaternary + start + 1)
    
    print("Calculating quaternary redundancy for strands up to {} bits for {}".format(strand_requirements.MAX_n_binary,
                                                                                      method))
    
    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary + start)):
        redundancy[n] = temp_strand_redundancy_array[n, common.COL_QUATERNARY_REDUNDANCY]

    redundancy = pd.DataFrame(redundancy)
    redundancy.index.name = "Original Quaternary Strand Length"
    redundancy.columns = [method + " quaternary redundancy"]
    redundancy.to_csv("EncodingAlgorithms/Results/" + method + "- Quaternary Redundancy.csv")
    return redundancy


def compute_binary_redundancy(strand_redundancy_df, requirements_name):
    binary_redundancy_arr = np.zeros(strand_requirements.MAX_n_binary + 1)

    start = 1
    optimal_n_quaternary = 0

    print("Calculating binary redundancy for strands up to {} bits for {}".format(strand_requirements.MAX_n_binary,
                                                                           requirements_name))
    for n in tqdm(range(start, strand_requirements.MAX_n_binary+start)):
        # calculate the first length that may represent a binary string of length n:
        optimal_n_quaternary = strand_redundancy_df.loc[strand_redundancy_df['Original_length'] >= ceil(n/2.0),
                                                        'Final_length'].min()
        
        # when max_binary_length >= n, the length is sufficient to represent 2^n possible vectors.
        binary_redundancy = int(2*optimal_n_quaternary-n)

        binary_redundancy_arr[n] = binary_redundancy

    binary_redundancy_arr = pd.DataFrame(binary_redundancy_arr)
    binary_redundancy_arr.index.name = "Original Binary String Length"
    binary_redundancy_arr.columns = [method + " Redundancy in bits"]
    binary_redundancy_arr.to_csv("EncodingAlgorithms/Results/" + method + "- Binary Redundancy.csv")
    return binary_redundancy_arr


if __name__ == "__main__":
    log2_count_array = None
    quaternary_redundancy_array = None
    binary_redundancy_array = None
    encoding_method = ["3-RLL and Knuth", "3-RLL and Knuth with Skips", "3-RLL and Knuth by blocks"]
    
    for method in encoding_method:
        if (method == "3-RLL and Knuth"):
            temp_strand_redundancy_array = RLL3andKnuth.calc_encoding_redundancy()
        elif (method == "3-RLL and Knuth with Skips"):
            temp_strand_redundancy_array = RLL3andKnuthwithSkips.calc_encoding_redundancy()
        elif (method == "3-RLL and Knuth by blocks"):
            temp_strand_redundancy_array = RLL3andKnuthbyBlocks.calc_encoding_redundancy()
        else:
            raise
        
        temp_strand_redundancy_df = pd.DataFrame(temp_strand_redundancy_array)
        temp_strand_redundancy_df.columns = ["Quaternary_redundancy", "Final_length"]
        temp_strand_redundancy_df = temp_strand_redundancy_df.rename_axis('Original_length').reset_index()
        temp_log2_count_array = compute_log2_count(temp_strand_redundancy_df, method)
        temp_quaternary_redundancy_array = compute_quaternary_redundancy(temp_strand_redundancy_array, method)
        temp_binary_redundancy_array = compute_binary_redundancy(temp_strand_redundancy_df, method)

        if (log2_count_array is None):
            log2_count_array = temp_log2_count_array
            quaternary_redundancy_array = temp_quaternary_redundancy_array
            binary_redundancy_array = temp_binary_redundancy_array
        else:
            log2_count_array = np.append(log2_count_array, temp_log2_count_array, axis=1)
            quaternary_redundancy_array = np.append(quaternary_redundancy_array, temp_quaternary_redundancy_array,
                                                    axis=1)
            binary_redundancy_array = np.append(binary_redundancy_array, temp_binary_redundancy_array, axis=1)


    log2_count_array = pd.DataFrame(log2_count_array, columns=encoding_method)
    log2_count_array.index.name = "Final Quaternary Strand Length"
    log2_count_array.to_csv("EncodingAlgorithms/Results/General Log2 Count.csv")

    quaternary_redundancy_array = pd.DataFrame(quaternary_redundancy_array, columns=encoding_method)
    quaternary_redundancy_array.index.name = "Original Quaternary Strand Length"
    quaternary_redundancy_array.to_csv("EncodingAlgorithms/Results/General Quaternary Redundancy.csv")

    binary_redundancy_array = pd.DataFrame(binary_redundancy_array, columns=encoding_method)
    binary_redundancy_array.index.name = "Original Binary String Length"
    binary_redundancy_array.to_csv("EncodingAlgorithms/Results/General Binary Redundancy.csv")
