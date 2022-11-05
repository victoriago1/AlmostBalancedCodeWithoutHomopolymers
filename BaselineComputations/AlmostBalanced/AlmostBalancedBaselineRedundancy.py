from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd
from pyrsistent import m
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
import strand_requirements

"""
Calculates the log2 of the number of almost balanced quaternary strands for every strand length n, based on
parameters from file strand_requirements.
This is done by summing combinatorically the number of strands.
In the created CSV, the columns are the value of n and its Log2 Count.
"""

def calc_strands_count():
    """
    Returns an array, for every n (length of quaternary strand) the value is the number of possible strands of
    length n that satisfy the constraint. """

    start = 2
    arr = np.empty(strand_requirements.MAX_n_quaternary + 1)
    arr[0] = 0
    arr[1] = 2

    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary + 1)):
        min_w, max_w = strand_requirements.min_max_weight(n)
        sum = 0
        for i in range(min_w, max_w + 1, 1):
            sum += comb(n, i)

        # sum * 2^n = number of possible strands of length n
        # by log rules: log2(sum * 2^n) = log2(sum) + log2(2^n) = log2(sum) + n
        arr[n] = log2(sum) + n
    return arr


if __name__ == "__main__":
    arr = calc_strands_count()
    df = pd.DataFrame(arr, columns=['Log2 Count'])
    df.index.name = "Final Quaternary String Length"
    df.to_csv("BaselineComputations/Results/Almost Balanced Baseline Redundancy.csv")
