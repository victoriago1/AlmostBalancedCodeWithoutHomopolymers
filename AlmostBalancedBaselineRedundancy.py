from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd

if __name__ == "__main__":
    N = 1000
    start = 2
    arr = np.empty((N,2))
    for n in range(start, N+start):
        arr_n_index = n - start
        print(n)
        sum = np.sum([comb(n, i) for i in range(ceil(0.45*n), floor(0.55*n) + 1, 1)])
        log_sum = floor(log2(sum)) if sum>0 else floor(log2(comb(n, n//2)))  # in case of odd n, by definition of exactly balanced (Knuth)
        arr[arr_n_index][0] = n
        arr[arr_n_index][1] = n - log_sum
    
    pd.DataFrame(arr).to_csv("./redundancy_of_almost_balances_constraint.csv")