import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
import strand_requirements

def min_encoded_length_by_method(x1, x2, x3, x4, x5):

    ''' calculates min length of letters replacing original zero runs and subtracts the original zeros length '''
    methods = {'A' : (2*(x1+x2+x3) + 3*(x4+x5)),
               'B' : (x1 + 3*(x2+x3) + 4*(x4+x5)),
               'C' : (x1 + 2*x2 + 6*(x3+x4+x5) + 1)}
    
    return min(methods, key=methods.get), (min(methods.values()) - (x1+2*x2+3*x3+4*x4+5*x5))
    # TODO return vector of mins

def compute_redundancy():
    methods = {'A' : 0,
               'B' : 0,
               'C' : 0}

    start = 80
    max_n = strand_requirements.MAX_n_quaternary
    results_per_n = np.empty((max_n+1) - start)

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
    
        results_per_n[n-start] = max_redundancy_per_length

    df = pd.DataFrame(results_per_n, columns=['max redundancy in quaternary letters'])
    df.to_csv("EncodingAlgorithms/Results/max_quaternary_redundancy_for_encoding_versions_idea.csv")

    df = pd.DataFrame.from_dict(methods, orient='index', columns=['counter per method'])
    df.to_csv("EncodingAlgorithms/Results/method_counter_for_encoding_versions_idea.csv")

    return results_per_n

if __name__ == "__main__":
    compute_redundancy()