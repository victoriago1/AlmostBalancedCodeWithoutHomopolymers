from math import log2
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
import strand_requirements

"""
Calculates the number of possible strands for the constraint of no homopolymer runs of length > 3,
without the constraint of balance or almost balance.
This is done by the recursive function: F(n) = 3*(F(n-1) + F(n-2) + F(n-3)), where F(n) = 4^n for n=1,2,3.
This file doesn't use the "m" weight variable from BaselineComputations\strand_requirements.py, because
the number of strands is calculated by a custom recursice function.
"""

def calc_strands_count():
    """
    Returns an array, for every n (length of quaternary strand) the value is the number of possible strands of
    that length that satisfy the constraint. """
    
    start = 1
    arr = np.zeros(strand_requirements.MAX_n_quaternary + 1)
    arr[0] = 0
    prev1 = prev2 = prev3 = 0

    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary+start)):
        num_of_strands = 4 ** n if n<=3 else 3 * (prev1 + prev2 + prev3)
        arr[n] = log2(num_of_strands)

        prev3 = prev2
        prev2 = prev1
        prev1 = num_of_strands
    return arr

if __name__ == "__main__":
    arr = calc_strands_count()
    df = pd.DataFrame(arr, columns=['Log2 Count'])
    df.index.name = "Final Quaternary String Length"
    df.to_csv("BaselineComputations/Results/3-RLL Baseline Redundancy.csv")
