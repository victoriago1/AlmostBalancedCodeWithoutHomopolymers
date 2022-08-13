import numpy as np
import pandas as pd
import sys
sys.path.append('./BaselineComputations/')
from strand_requirements import  m, MAX_n, MAX_w

class PowerSeries:
    def __init__(self, rows=-1, cols=-1):
        if (rows > 0):
            self.is_compact = False
            self.data = np.zeros((rows, cols))
        else:
            self.is_compact = None
            self.data = None
    
    def init_T(self): # This initiates the power series according to the function T defined in the article
        self.is_compact = False
        self.data = np.ones((m+1,1))
        self.data[0,0] = 0

    def init_T1(self):
        self.is_compact = True
        self.data = np.ones(m+1)
        self.data[0] = 0


    def __mul__(self, other):
        ps_self = np.diag(self.data) if self.is_compact else self.data
        ps_other = np.diag(other.data) if other.is_compact else other.data

        new_rows = min(ps_self.shape[0] + ps_other.shape[0], MAX_n + 1)
        new_cols = min(ps_self.shape[1] + ps_other.shape[1], MAX_w + 1)
        result = PowerSeries(new_rows, new_cols)

        for self_i in range(min(ps_self.shape[0], new_rows)):
            for self_j in range(min(ps_self.shape[1], new_cols)):
                
                self_coeff = ps_self[self_i, self_j]

                for other_i in range(min(ps_other.shape[0], new_rows-self_i)):
                    for other_j in range(min(ps_other.shape[1], new_cols-self_j)):
                        result.data[self_i+other_i, self_j+other_j] += self_coeff*ps_other[other_i, other_j]
    
        return result


    def add_matrix_to_self(self, ps_matrix):

        # note: self.is_compact must be false
        # note: ps_other shape should be in this point smaller or equal to self

        for matrix_i in range(ps_matrix.shape[0]):
            for matrix_j in range(ps_matrix.shape[1]):
                self.data[matrix_i, matrix_j] += ps_matrix[matrix_i, matrix_j]
        
        return self
    
    def __add__(self, other):
        ps_self = np.diag(self.data) if self.is_compact else self.data
        ps_other = np.diag(other.data) if other.is_compact else other.data
        
        new_rows = max(ps_self.shape[0], ps_other.shape[0])
        new_cols = max(ps_self.shape[1], ps_other.shape[1])
        result = PowerSeries(new_rows, new_cols)

        result.add_matrix_to_self(ps_self)

        return result.add_matrix_to_self(ps_other)


    def __str__(self):
        str = ""
        if self.is_matrix:
            ps_self = self.data if (len(self.data.shape) > 1) else np.diag(self.data)
            # matrix, x and y
            # TODO: fix efficency 
            for (i,j), val in np.ndenumerate(ps_self):
                if val:
                    if (j==0):
                        str += '{}*(x^{}) + '.format(val, i)
                    else:
                        str += '{}*(x^{})*(y^{}) + '.format(val, i, j)
        else:
            # vector, only x
            for i in range(len(self.data)):
                if self.data[i]:
                    str += '{}*(x^{}) + '.format(self.data[i], i)
        return str[:-3] if len(str) else "0"