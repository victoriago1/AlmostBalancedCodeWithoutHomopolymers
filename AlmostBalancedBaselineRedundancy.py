from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm

"""
Calculates the number of bits for minimum and optimal redundancy, in case of almost balanced code,
without the constraint of GC/AT content.
This is done by counting the number of the strands that hold this condition, and calculating the
number of binary bits that can express the same number of strands, and returning the difference between the two.
In the created CSV, the columns are the value of n and its redundancy.
"""

if __name__ == "__main__":
    N = 1000
    start = 2
    arr = np.empty((N,2))
    for n in tqdm(range(start, N+start)):
        arr_n_index = n - start
        sum = np.sum([comb(n, i) for i in range(ceil(0.45*n), floor(0.55*n) + 1, 1)])
        # weight is exactly n//2 in case of odd n and small values (and then sum == 0), by definition of exactly balanced (Knuth)
        log_sum = floor(log2(sum)) if sum>0 else floor(log2(comb(n, n//2))) 
        arr[arr_n_index][0] = n
        arr[arr_n_index][1] = n - log_sum
    
    pd.DataFrame(arr).to_csv("./redundancy_of_almost_balances_constraint.csv")