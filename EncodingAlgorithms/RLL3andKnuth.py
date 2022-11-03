from math import log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
sys.path.append('BaselineComputations/RLL/')
sys.path.append('EncodingAlgorithms/')
import strand_requirements
import RLL3BaselineRedundancy
import common

"""
Calculates the number of possible strands for the algorithm of 3-RLL & Knuth.
"""

def calc_encoding_redundancy():
    '''
    Returns an array, for every n (length of quaternary strand) the value is the length of the strand that satisfy the
    constraint after encoding. '''
    T2_log_count_array = RLL3BaselineRedundancy.calc_strands_count()
    start = 1
    arr = np.zeros((strand_requirements.MAX_n_quaternary + 1, 2))

    arr[0, common.COL_QUATERNARY_REDUNDANCY] = 0
    arr[0, common.COL_FINAL_LENGTH] = 0
    index_size = 0
    
    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary+start)):
        n_float = 1.0*n
        RLL_redundancy = 2*ceil(n_float/16) - 1
        n_length_after_RLL = n + RLL_redundancy

        while T2_log_count_array[index_size] < log2(n_length_after_RLL):
            index_size += 1
       
        Knuth_redundancy = index_size + 2
        arr[n, common.COL_QUATERNARY_REDUNDANCY] = RLL_redundancy + Knuth_redundancy
        arr[n, common.COL_FINAL_LENGTH] = n + arr[n, common.COL_QUATERNARY_REDUNDANCY]

    return arr

if __name__ == "__main__":
    arr = calc_encoding_redundancy()
    df = pd.DataFrame(arr, columns=['Quaternary Redundancy', 'Final Quaternary String Length'])
    df.index.name = "Original Quaternary String Length"
    df.to_csv("EncodingAlgorithms/Results/3-RLL and Knuth encoding.csv")
