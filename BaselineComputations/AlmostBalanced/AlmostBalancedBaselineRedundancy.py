from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd
from pyrsistent import m
from tqdm import tqdm
import sys
sys.path.append('./BaselineComputations/')
import strand_requirements

"""
Calculates the number of bits for minimum and optimal redundancy, in case of almost balanced code,
without the constraint of GC/AT content.
This is done by
(1) counting the number of the 4-ary strands of length n that hold this condition
(2) calculating the number of binary bits that can express the same number of strands
(3) returning the difference between 2n (the binary bits that are equal to n 4-ary strands) and the result of (2)
In practice we used equations simplification with log rules.
In the created CSV, the columns are the value of n and its redundancy.
"""

def calc_strands_count():
    '''
    Returns an array, for every n (length of quaternary strand) the value is the number of possible strands of
    that length that satisfy the constraint. '''
    
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
        arr[n] = log2(sum) + n
    return arr


if __name__ == "__main__":
    arr = calc_strands_count()
    df = pd.DataFrame(arr, columns=['max length of binary vectors'])
    df.to_csv("./BaselineComputations/Results/max_binary_vectors_length_for_strictly_balances_constraint.csv")
