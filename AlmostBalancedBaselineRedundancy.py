import math
import numpy as np

if __name__ == "__main__":
    arr = np.array(1000)
    for n in range(1000):
        sum = np.sum([math.comb(n, i) for i in range(math,ceil)])