from math import comb, log2, ceil, floor
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

COL_INDEX_SIZE = 0
COL_DELTA = 1

def compute_index_size_and_delta (max_n):
    T2_log_count_array = RLL3BaselineRedundancy.calc_strands_count()
    imbalance = 0.01*strand_requirements.IMBALANCE_PERCENTAGE
    arr = np.empty((max_n + 1, 2))
    arr[0, COL_INDEX_SIZE] =  arr[0, COL_DELTA] = 0
    arr[1, COL_INDEX_SIZE] = 1
    arr[1, COL_DELTA] = 2
    index_size = 0
    for n in range(2, max_n + 1):
        n_float= 1.0*n

        while T2_log_count_array[index_size] < log2(n_float):
                index_size += 1
        
        arr[n, COL_INDEX_SIZE] = index_size
        arr[n, COL_DELTA] = 1
        for delta in range(ceil(0.1*(n+index_size+2))):
            while T2_log_count_array[index_size-1] >= log2(n_float/(delta+1)):
                index_size -= 1
            if ((imbalance*(n_float+index_size+2)))>=ceil(0.5*(delta+1+index_size+2)):
                arr[n, COL_INDEX_SIZE] = index_size
                arr[n, COL_DELTA] = delta+1
    
    df = pd.DataFrame(arr, columns=['Index Size', 'Optimal Delta'])
    df.index.name = "Quaternary String Length"
    df.to_csv("EncodingAlgorithms/Results/Index and optimal delta per length.csv")
    return arr

def calc_encoding_redundancy():
    '''
    Returns an array, for every n (length of quaternary strand) the value is the length of the strand that satisfy the
    constraint after encoding. '''

    max_length_for_delta_calc = ceil(1.25*strand_requirements.MAX_n_quaternary + 1)
    delta_and_index = compute_index_size_and_delta(max_n=max_length_for_delta_calc)
    start = 1
    arr = np.zeros((strand_requirements.MAX_n_quaternary + 1, 2))

    arr[0, common.COL_QUATERNARY_REDUNDANCY] = 0
    arr[0, common.COL_FINAL_LENGTH] = 0
    
    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary+start)):
        n_float = 1.0*n
        RLL_redundancy = 2*ceil(n_float/16) - 1
        n_length_after_RLL = n + RLL_redundancy
        index_size = delta_and_index[n_length_after_RLL, COL_INDEX_SIZE]
        Knuth_redundancy = index_size + 2
        arr[n, common.COL_QUATERNARY_REDUNDANCY] = RLL_redundancy + Knuth_redundancy
        arr[n, common.COL_FINAL_LENGTH] = n + arr[n, common.COL_QUATERNARY_REDUNDANCY]
        if  arr[n, common.COL_FINAL_LENGTH] > max_length_for_delta_calc:
            break

    return arr

if __name__ == "__main__":
    arr = calc_encoding_redundancy()
    df = pd.DataFrame(arr, columns=['Quaternary Redundancy', 'Final Quaternary String Length'])
    df.index.name = "Original Quaternary String Length"
    df.to_csv("EncodingAlgorithms/Results/3-RLL and Knuth with skips encoding.csv")
