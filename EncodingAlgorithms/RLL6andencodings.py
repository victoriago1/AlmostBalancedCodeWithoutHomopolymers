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
Calculations for the third proposed solution algorithm, with 6-RLL & methods A, B, C.
Currently, Computed only for values above 80.
"""

def min_encoded_length_by_method(x1, x2, x3, x4, x5):
    """ Calculates min of length of letters replacing original zero runs and subtracts the original zeros length.
        This equals to the redundancy for the given zero runs frequencies in a strand. """

    methods = {'A' : (2*(x1+x2+x3) + 3*(x4+x5)),
               'B' : (x1 + 3*(x2+x3) + 4*(x4+x5)),
               'C' : (x1 + 2*x2 + 6*(x3+x4+x5) + 1)}

    # Note that when two methods have the same redundancy, the first one will be chosen.
    # This affects the choosing frequency of each method.
    return min(methods, key=methods.get), (min(methods.values()) - (x1+2*x2+3*x3+4*x4+5*x5))


def calc_encoding_redundancy():
    """
    Returns an array, for every n (length of quaternary strand) the value is the worst-case redundancy for that length
    with the algorithm, among all three encoding methods. """

    methods = {'A' : 0,
               'B' : 0,
               'C' : 0}

    max_n = strand_requirements.MAX_n_quaternary
    delta_and_index = common.compute_index_size_and_delta(max_n=(max_n+ceil(0.5*max_n)))
    start = 80
    results_per_n = np.empty(((max_n+1), 2))

    for n in tqdm(range(start, max_n+1)):
        max_redundancy_per_length = 0
        for t1 in range(0,(n+2),2):
            x1 = t1/2
            for t2 in range(0,(n+2)-t1,3):
                x2 = t2/3
                for t3 in range(0,(n+2)-t1-t2,4):
                    x3 = t3/4
                    for t4 in range(0,(n+2)-t1-t2-t3,5):
                        x4 = t4/5
                        
                        t5 = ((n+1)-t1-t2-t3-t4)
                        if (t5%6 != 0):
                            continue
                        
                        x5 = t5/6
                        method, min_redundancy_per_strand = min_encoded_length_by_method(x1, x2, x3, x4, x5)

                        methods[method] += 1
                        
                        if (min_redundancy_per_strand > max_redundancy_per_length):
                            max_redundancy_per_length = min_redundancy_per_strand
    
        RLL_redundancy = ceil(max_redundancy_per_length)
        n_length_after_RLL = n + RLL_redundancy
        index_size = delta_and_index[n_length_after_RLL, common.COL_INDEX_SIZE]
        Knuth_redundancy = index_size + 2
        results_per_n[n, common.COL_QUATERNARY_REDUNDANCY] = RLL_redundancy + Knuth_redundancy
        results_per_n[n, common.COL_FINAL_LENGTH] = n + results_per_n[n, common.COL_QUATERNARY_REDUNDANCY]

    df = pd.DataFrame.from_dict(methods, orient='index', columns=['counter per method']) # note the comment above
    df.to_csv("EncodingAlgorithms/Results/6-RLL with encodings and Knuth - counter for zero runs encoding methods.csv") 

    return results_per_n


if __name__ == "__main__":
    arr = calc_encoding_redundancy()
    df = pd.DataFrame(arr, columns=['Quaternary Redundancy', 'Final Quaternary String Length'])
    df.index.name = "Original Quaternary String Length"
    df.to_csv("EncodingAlgorithms/Results/6-RLL with encodings and Knuth encoding.csv")
