from math import log2, floor
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('BaselineComputations/')
import strand_requirements
from PowerSeries import PowerSeries


def _create_D1():
    D1 = np.full((4, 4), None)

    for (i,j), val in np.ndenumerate(D1):
        if (i==j):
            D1[i,j] = PowerSeries(1, 1)
        
        elif (i<=1):
            D1[i,j] = PowerSeries()
            D1[i,j].init_T()
        
        else:
            D1[i,j] = PowerSeries()
            D1[i,j].init_T1()

    return D1    

def _print_D(D):
    for (i,j), val in np.ndenumerate(D):
        if(j == 0):
            print("*******************************************")
        print("i,j = {},{}: {}".format(i, j, val))

def calc_strands_count():
    D1 = _create_D1()

    N4 = PowerSeries(strand_requirements.MAX_n_quaternary + 1, strand_requirements.MAX_w + 1)

    print("Calculating N4 function up till {}".format(strand_requirements.MAX_n_quaternary))
    for i in tqdm(range(strand_requirements.MAX_n_quaternary)):
        if(i==0):
            Dk = _create_D1()
        else:
            Dk = Dk@D1

        for (i,j), val in np.ndenumerate(Dk):
            N4 += val
    
    num_of_strands = np.zeros((strand_requirements.MAX_n_quaternary + 1, 2), dtype=np.float128)
    COL_NUM_OF_STRANDS = 0
    COL_MAX_BINARY_LENGTH = 1
    num_of_strands[0] = 0

    for n in range(1, strand_requirements.MAX_n_quaternary + 1):
        min_w, max_w = strand_requirements.min_max_weight(n)
        for w in range(min_w, max_w+1):
            num_of_strands[n, COL_NUM_OF_STRANDS] += (1/3)*N4.data[n,w]

        num_of_strands[n, COL_MAX_BINARY_LENGTH] = log2(num_of_strands[n, COL_NUM_OF_STRANDS], dtype=np.float128)

    pd.DataFrame(N4.data).to_csv("BaselineComputations/Results/Almost Balanced And RLL- N4_PowerSeries.csv")

    strands_count = num_of_strands   
    num_of_strands = pd.DataFrame(num_of_strands)
    num_of_strands.index.name = "Final Quaternary String Length"
    num_of_strands.columns = ["Number of possible strands", "Log2 Count"]
    num_of_strands.to_csv("BaselineComputations/Results/Almost Balanced And RLL Baseline Redundancy (N4 calc).csv")

    return strands_count[:,COL_MAX_BINARY_LENGTH]

if __name__ == "__main__":
    
    calc_strands_count()