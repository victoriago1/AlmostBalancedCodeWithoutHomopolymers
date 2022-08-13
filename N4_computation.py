from math import comb, fabs, log2, ceil, floor
import numpy as np
import pandas as pd
import common
from PowerSeries import PowerSeries
from tqdm import tqdm


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

    N4 = PowerSeries(common.MAX_n + 1, common.MAX_w + 1)

    for i in tqdm(range(common.MAX_n)):
        if(i==0):
            Dk = _create_D1()
        else:
            Dk = Dk@D1

        for (i,j), val in np.ndenumerate(Dk):
            N4 += val
    
    num_of_strands = [0]
    for n in range(1, common.MAX_n + 1):
        num_of_strands.append(0)
        min_w = ceil(0.45*n)
        max_w = floor(0.55*n)
        max_w = max_w + 1 if (max_w<min_w) else max_w
        for w in range(min_w, max_w+1):
            num_of_strands[n]+=(1/3)*N4.data[n,w]

    print(num_of_strands)
    pd.DataFrame(N4.data).to_csv("./N4_calc.csv")
    pd.DataFrame(num_of_strands).to_csv("./number_of_strands_by_length_using_N4.csv")