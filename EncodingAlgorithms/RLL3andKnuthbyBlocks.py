from math import ceil
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
sys.path.append('EncodingAlgorithms/')
import strand_requirements
import common

"""
Calculates the number of possible strands for the algorithm of 3-RLL & Knuth by blocks (our leading idea).
"""

def calc_encoding_redundancy():
    '''
    Returns an array, for every n (length of quaternary strand) the value is the length of the strand that satisfy the
    constraint after encoding. '''
    start = 1
    arr = np.zeros((strand_requirements.MAX_n_quaternary + 1, 2))

    arr[0, common.COL_QUATERNARY_REDUNDANCY] = 0
    arr[0, common.COL_FINAL_LENGTH] = 0
    
    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary+start)):
        n_float = 1.0*n
        RLL_redundancy = ceil(n_float/16)
        Knuth_redundancy = 2 if n_float >= 186 else 3 
        arr[n, common.COL_QUATERNARY_REDUNDANCY] = RLL_redundancy + Knuth_redundancy
        arr[n, common.COL_FINAL_LENGTH] = n + arr[n, common.COL_QUATERNARY_REDUNDANCY]

    return arr

if __name__ == "__main__":
    arr = calc_encoding_redundancy()
    df = pd.DataFrame(arr, columns=['Quaternary Redundancy', 'Final Quaternary String Length'])
    df.index.name = "Original Quaternary String Length"
    df.to_csv("EncodingAlgorithms/Results/3-RLL and Knuth by blocks (our idea) encoding.csv")