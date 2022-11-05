from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
sys.path.append('EncodingAlgorithms/')
import strand_requirements
import common


"""
Calculations for the first proposed solution algorithm, with 3-RLL & Knuth with skips.
"""

def calc_encoding_redundancy():
    """
    Returns an array, for every n (length of quaternary strand) the value is the length of the strand that satisfies the
    constraint after encoding, and the total encoding redundancy for that length. """

    max_length_for_delta_calc = ceil(1.25*strand_requirements.MAX_n_quaternary + 1)
    delta_and_index = common.compute_index_size_and_delta(max_n=max_length_for_delta_calc)
    start = 1
    arr = np.zeros((strand_requirements.MAX_n_quaternary + 1, 2))

    arr[0, common.COL_QUATERNARY_REDUNDANCY] = 0
    arr[0, common.COL_FINAL_LENGTH] = 0

    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary+start)):
        n_float = 1.0*n
        RLL_redundancy = 2*ceil(n_float/16) - 1
        n_length_after_RLL = n + RLL_redundancy
        index_size = delta_and_index[n_length_after_RLL, common.COL_INDEX_SIZE]
        Knuth_redundancy = index_size + 2
        arr[n, common.COL_QUATERNARY_REDUNDANCY] = RLL_redundancy + Knuth_redundancy
        arr[n, common.COL_FINAL_LENGTH] = n + arr[n, common.COL_QUATERNARY_REDUNDANCY]
        if arr[n, common.COL_FINAL_LENGTH] > max_length_for_delta_calc:
            break

    return arr


if __name__ == "__main__":
    arr = calc_encoding_redundancy()
    df = pd.DataFrame(arr, columns=['Quaternary Redundancy', 'Final Quaternary String Length'])
    df.index.name = "Original Quaternary String Length"
    df.to_csv("EncodingAlgorithms/Results/3-RLL and Knuth with skips encoding.csv")
