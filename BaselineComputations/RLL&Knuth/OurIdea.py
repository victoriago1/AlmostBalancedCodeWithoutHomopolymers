from math import comb, log2, ceil, floor
import numpy as np
import pandas as pd
import sys
sys.path.append('./BaselineComputations/')
import strand_requirements

def compute_deltas (max_n, max_delta):
    arr = np.empty(max_n)
    arr[0] = 0
    for n in range(1,max_n):
        n_float= 1.0*n
        arr[n] = 0
        for delta in range(max_delta):
            index_size = ceil((log2(n_float/(delta+1)))/2)
            if ((0.05*(n_float+index_size)))>=ceil(0.5*(delta+1+index_size)):
                arr[n] = delta+1
            # else:
            #     break
    
    pd.DataFrame(arr).to_csv("BaselineComputations/Results/Max delta per length.csv")
    return arr

if __name__ == "__main__":
    # compute_deltas(max_n=strand_requirements.MAX_n, max_delta=ceil(0.025*strand_requirements.MAX_n))
    compute_deltas(max_n=1000, max_delta=ceil(0.1*1000))

