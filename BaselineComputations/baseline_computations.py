from math import log2, ceil, floor
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
import strand_requirements

sys.path.append('BaselineComputations/AlmostBalanced/')
sys.path.append('BaselineComputations/AlmostBalancedAndRLL/')
sys.path.append('BaselineComputations/RLL/')
import AlmostBalancedBaselineRedundancy
import N4_computation
import RLL3BaselineRedundancy

"""
A script for running all optimal bounds calculations and saving them in seperate CSV files.
"""

def compute_log2_count(strand_count_array, method):
    log2_count = strand_count_array
    log2_count = pd.DataFrame(log2_count)
    log2_count.index.name = "Final Quaternary Strand Length"
    log2_count.columns = [method + " log2 of possible strands"]
    log2_count.to_csv("BaselineComputations/Results/" + method + "- Log2 Count.csv")
    return log2_count


def compute_quaternary_redundancy(strand_count_array, method):
    start = 1
    redundancy = np.zeros(strand_requirements.MAX_n_quaternary + start)
    optimal_final_n_quaternary = max_origin_quaternary = 0
    print("Calculating quaternary redundancy for strands up to {} bits for {}".format(strand_requirements.MAX_n_binary,
                                                                                      method))

    for n in tqdm(range(start, strand_requirements.MAX_n_quaternary + start)):
        # calculate the first length that may represent a binary string of length n:
        while max_origin_quaternary < n:
            optimal_final_n_quaternary += 1
            max_origin_quaternary = floor(strand_count_array[optimal_final_n_quaternary]/2)
        
        redundancy[n] = optimal_final_n_quaternary - n

        if optimal_final_n_quaternary >= strand_requirements.MAX_n_quaternary:
            break

    redundancy = pd.DataFrame(redundancy)
    redundancy.index.name = "Original Quaternary Strand Length"
    redundancy.columns = [method + " quaternary redundancy"]
    redundancy.to_csv("BaselineComputations/Results/" + method + "- Quaternary Redundancy.csv")
    return redundancy


def compute_binary_redundancy(max_n_binary_per_n_quaternary, requirements_name):
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
    redundancy.index.name = "Original Binary String Length"
    redundancy.columns = [requirements_name + " redundancy in bits"]
    redundancy.to_csv("BaselineComputations/Results/" + requirements_name + "- Binary Redundancy.csv")
    return redundancy


if __name__ == "__main__":
    log2_count_array = None
    quaternary_redundancy_array = None
    binary_redundancy_array = None
    requirements = ["Almost Balanced", "3-RLL", "Almost Balanced And RLL"]
    
    for requirement in requirements:
        if (requirement == "Almost Balanced"):
            temp_strand_count_array = AlmostBalancedBaselineRedundancy.calc_strands_count()
        elif (requirement == "3-RLL"):
            temp_strand_count_array = RLL3BaselineRedundancy.calc_strands_count()
        elif (requirement == "Almost Balanced And RLL"):
            temp_strand_count_array = N4_computation.calc_strands_count()
        else:
            raise
        
        temp_log2_count_array = compute_log2_count(temp_strand_count_array, requirement)
        temp_quaternary_redundancy_array = compute_quaternary_redundancy(temp_strand_count_array, requirement)
        temp_binary_redundancy_array = compute_binary_redundancy(temp_strand_count_array, requirement)

        if (log2_count_array is None):
            log2_count_array = temp_log2_count_array
            quaternary_redundancy_array = temp_quaternary_redundancy_array
            binary_redundancy_array = temp_binary_redundancy_array
        else:
            log2_count_array = np.append(log2_count_array, temp_log2_count_array, axis=1)
            quaternary_redundancy_array = np.append(quaternary_redundancy_array, temp_quaternary_redundancy_array,
                                                    axis=1)
            binary_redundancy_array = np.append(binary_redundancy_array, temp_binary_redundancy_array, axis=1)


    log2_count_array = pd.DataFrame(log2_count_array, columns=requirements)
    log2_count_array.index.name = "Final Quaternary Strand Length"
    log2_count_array.to_csv("BaselineComputations/Results/General Log2 Count.csv")

    quaternary_redundancy_array = pd.DataFrame(quaternary_redundancy_array, columns=requirements)
    quaternary_redundancy_array.index.name = "Original Quaternary Strand Length"
    quaternary_redundancy_array.to_csv("BaselineComputations/Results/General Quaternary Redundancy.csv")

    binary_redundancy_array = pd.DataFrame(binary_redundancy_array, columns=requirements)
    binary_redundancy_array.index.name = "Original Binary String Length"
    binary_redundancy_array.to_csv("BaselineComputations/Results/General Binary Redundancy.csv")
