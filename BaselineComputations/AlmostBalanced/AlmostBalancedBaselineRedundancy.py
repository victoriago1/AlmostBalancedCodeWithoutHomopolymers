from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm

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

if __name__ == "__main__":
    N = 1000
    start = 1
    arr = np.empty((N,2))

    COL_N = 0
    COL_REDUNDANCY = 1

    for n in tqdm(range(start, N+start)):
        arr_n_index = n - start
        sum = np.sum([comb(n, i) for i in range(ceil(0.45*n), floor(0.55*n) + 1, 1)])
        # weight is exactly n//2 in case of odd n and small values (and then sum == 0), by definition of exactly balanced (Knuth)
        log_sum = floor(log2(sum)) if sum>0 else floor(log2(comb(n, n//2))) 
        arr[arr_n_index][COL_N] = n
        arr[arr_n_index][COL_REDUNDANCY] = n - log_sum
    
    df = pd.DataFrame(arr, columns=['4-ary n', 'binary redundancy'])
    df.to_csv("./redundancy_of_almost_balances_constraint.csv")