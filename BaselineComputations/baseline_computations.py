from math import log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
import strand_requirements

sys.path.append('./BaselineComputations/AlmostBalanced/')
sys.path.append('./BaselineComputations/AlmostBalancedAndRLL/')
sys.path.append('./BaselineComputations/RLL/')
import AlmostBalancedBaselineRedundancy
import N4_computation
import NoHomopolymersBaselineRedundancy

def compute_redundancy(max_n_binary_per_n_quaternary, requirements_name):
    redundancy = np.zeros(strand_requirements.MAX_n_binary + 1)

    start = 1
    optimal_n_quaternary = ceil(start/2.0)
    max_binary_length = max_n_binary_per_n_quaternary[optimal_n_quaternary]
    

    print("Calculating redundancy for strands up to {} bits for {}".format(strand_requirements.MAX_n_binary,
                                                                           requirements_name))
    for n in tqdm(range(start, strand_requirements.MAX_n_binary+start)):
        # calculate the first length that may represent a binary string of length n:
        while max_binary_length < n:
            optimal_n_quaternary += 1
            max_binary_length = max_n_binary_per_n_quaternary[optimal_n_quaternary]
        
        # when max_binary_length >= n, the length is sufficient to represent 2^n possible vectors.
        binary_redundancy = int(2*optimal_n_quaternary - n)

        redundancy[n] = binary_redundancy

    redundancy = pd.DataFrame(redundancy)
    redundancy.index.name = "Binary String Length"
    redundancy.columns = [requirements_name + " Redundancy in bits"]
    redundancy.to_csv("BaselineComputations/Results/" + requirements_name + "-Redundancy.csv")
    return redundancy

if __name__ == "__main__":
    redundancy_array = None
    strand_count_array = None
    requirements = ["AlmostBalanced", "RLL", "AlmostBalancedAndRLL"]
    
    for requirement in requirements:
        if (requirement == "AlmostBalanced"):
            temp_strand_count_array = AlmostBalancedBaselineRedundancy.calc_strands_count()
            temp_redundancy_array = compute_redundancy(temp_strand_count_array, requirement)
        elif (requirement == "RLL"):
            temp_strand_count_array = NoHomopolymersBaselineRedundancy.calc_strands_count()
            temp_redundancy_array =  compute_redundancy(temp_strand_count_array, requirement)
        elif (requirement == "AlmostBalancedAndRLL"):
            temp_strand_count_array = N4_computation.calc_strands_count()
            temp_redundancy_array = compute_redundancy(temp_strand_count_array, requirement)
        else:
            raise
        
        if (redundancy_array is None):
            redundancy_array = temp_redundancy_array
            strand_count_array = temp_strand_count_array.reshape(-1,1)
        else:
            redundancy_array = np.append(redundancy_array, temp_redundancy_array, axis=1)
            strand_count_array = np.append(strand_count_array, temp_strand_count_array.reshape(-1,1), axis=1)
        
        strand_count_array_df = pd.DataFrame(temp_strand_count_array)
        strand_count_array_df.index.name = "Quaternary String Length"
        strand_count_array_df.columns = [requirement + " log2 of possible strands count"]
        strand_count_array_df.to_csv("BaselineComputations/Results/" + requirement + "- Log2 Count.csv")

    redundancy_array = pd.DataFrame(redundancy_array, columns=requirements)
    redundancy_array.index.name = "Binary String Length"
    redundancy_array.to_csv("BaselineComputations/Results/General Redundancy.csv")

    strand_count_array = pd.DataFrame(strand_count_array, columns=requirements)
    strand_count_array.index.name = "Quaternary String Length"
    strand_count_array.to_csv("BaselineComputations/Results/General log2 strand count.csv")