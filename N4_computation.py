from math import comb, fabs, log2, ceil, floor
import numpy as np
import pandas as pd



m = 3
MAX_n = 1000
MAX_w = floor(0.55*MAX_n)

class PowerSeries:
    def __init__(self, rows = -1, cols=-1):
        if (rows > 0):
            self.is_matrix = True
            self.data  = np.zeros((rows, cols))
    
    def init_T(self):
        self.is_matrix = False
        self.data  = np.ones(m)

    def init_T1(self):
        self.is_matrix = True
        self.data  = np.ones(m)
    
    def __mul__(self, other):
        ps_self = self.data if (not self.is_matrix or self.data.shape[1]>0) else np.diag(self.data)
        ps_other = other.data if (not other.is_matrix or other.data.shape[1]>0) else np.diag(other.data)

        if (self.is_matrix and other.is_matrix):
            new_rows = min(self.data.shape[0] + other.data.shape[0], MAX_n)
            new_cols = min(self.data.shape[1] + other.data.shape[1], MAX_w)

            result = PowerSeries(new_rows, new_cols)

            # #option 1:
            # for i in range(new_rows):
            #     for j in range(new_cols):
                
            #         for k in range(i):
            #             for l in range(j):
            #                 if (k>(self.data.shape[0]-1) or (i-k)>(other.data.shape[0]-1)):
            #                     continue
            #                 if (l>(self.data.shape[1]-1) or (j-l)>(other.data.shape[1]-1)):
            #                     continue
            #                 result[i,j] += self.data[k,l]*other.data[i-k,j-l]
            
            #option 2:
            for self_i in range(min(ps_self.shape[0], new_rows)):
                for self_j in range(min(ps_self.shape[1], new_cols)):
                    
                    self_coeff = ps_self[self_i, self_j]

                    for other_i in range(min(ps_other.shape[0], new_rows-self_i)):
                        for other_j in range(min(ps_other.shape[1], new_cols-self_j)):
                            result[self_i+other_i, self_j+other_j] = self_coeff*ps_other[other_i, other_j]
            
        elif(self.is_matrix):
            for self_i in range(min(ps_self.shape[0], new_rows)):
                for self_j in range(min(ps_self.shape[1], new_cols)):
                    
                    self_coeff = ps_self[self_i, self_j]

                    for other_i in range(min(ps_other.shape[0], new_rows-self_i)):
                        result[self_i+other_i, self_j] = self_coeff*ps_other[other_i]

        else:
            for other_i in range(min(ps_other.shape[0], new_rows)):
                for other_j in range(min(ps_other.shape[1], new_cols)):                    
                   
                    other_coeff = ps_other[other_i, other_j]

                    for self_i in range(min(ps_self.shape[0], new_rows-other_i)):
                        result[self_i+other_i, other_j] = other_coeff*ps_self[self_i]
    
        return result
                


    

                
        

if __name__ == "__main__":