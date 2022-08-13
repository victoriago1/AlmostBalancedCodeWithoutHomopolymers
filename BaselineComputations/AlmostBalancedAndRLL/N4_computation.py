import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
sys.path.append('./BaselineComputations/')
import strand_requirements
from PowerSeries import PowerSeries


def _create_D1():
    D1 = np.full((4, 4), None)

    for (i,j), val in np.ndenumerate(D1):
        if (i==j):
            D1[i,j] = PowerSeries(1, 1, False)
        
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

if __name__ == "__main__":
    
    D1 = _create_D1()

    N4 = PowerSeries(strand_requirements.MAX_n + 1, strand_requirements.MAX_w + 1)

    for i in tqdm(range(strand_requirements.MAX_n)):
        if(i==0):
            Dk = _create_D1()
        else:
            Dk = Dk@D1

        for (i,j), val in np.ndenumerate(Dk):
            N4 += val
    
    num_of_strands = [0]
    for n in range(1, strand_requirements.MAX_n + 1):
        num_of_strands.append(0)
        min_w, max_w = strand_requirements.min_man_w(n)
        for w in range(min_w, max_w+1):
            num_of_strands[n]+=(1/3)*N4.data[n,w]

    print(num_of_strands)
    pd.DataFrame(N4.data).to_csv("BaselineComputations/Results/AlmostBalancedAndRLL-N4_PowerSeries.csv")
    num_of_strands = pd.DataFrame(num_of_strands)
    num_of_strands.index.name = "Strand length"
    num_of_strands.columns = ["Number of possible strands"]
    num_of_strands.to_csv("BaselineComputations/Results/AlmostBalancedAndRLL-number_of_strands_by_length_using_N4.csv")