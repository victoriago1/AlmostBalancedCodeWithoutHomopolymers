from math import log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm

"""
Calculates the number of bits for minimum and optimal redundancy, for the constraint of no homopolymers of length > 3,
without the constraint of balance or almost balance.
This is done by the recursice function: F(n) = 3*(F(n-1) + F(n-2) + F(n-3)), and F(n) = 4^n for n=1,2,3.
"""


if __name__ == "__main__":
    N = 1000
    start = 1
    arr = np.empty((N,2))
    prev1 = prev2 = prev3 = 0

    COL_N = 0
    COL_REDUNDANCY = 1

    for n in tqdm(range(start, N+start)):
        arr_n_index = n - start
        num_of_strands = 4 ** n if n<=3 else 3 * (prev1 + prev2 + prev3)
        max_binary_length = floor(log2(num_of_strands))
        binary_redundancy = 2 * n - max_binary_length

        arr[arr_n_index][COL_N] = n
        arr[arr_n_index][COL_REDUNDANCY] = binary_redundancy

        prev3 = prev2
        prev2 = prev1
        prev1 = num_of_strands
    
    df = pd.DataFrame(arr, columns=['4-ary n', 'binary redundancy'])
    df.to_csv("./redundancy_of_no_homopolymers_constraint.csv")