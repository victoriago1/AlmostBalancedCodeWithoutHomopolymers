from math import log2, ceil
import numpy as np
import pandas as pd
import sys
sys.path.append('BaselineComputations/')
sys.path.append('BaselineComputations/RLL/')
import strand_requirements
import RLL3BaselineRedundancy

COL_QUATERNARY_REDUNDANCY = 0
COL_FINAL_LENGTH = 1

COL_INDEX_SIZE = 0
COL_DELTA = 1

def compute_index_size_and_delta (max_n):
    """
    Calculates for every encoded strand length the optimal interval size (delta) and index encoding size.
    """
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