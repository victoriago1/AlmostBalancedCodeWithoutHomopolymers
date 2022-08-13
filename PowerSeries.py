import numpy as np
import pandas as pd
from common import  m, MAX_n, MAX_w

class PowerSeries:
    def __init__(self, rows=-1, cols=-1, is_matrix=True):
        if (rows > 0):
            self.is_matrix = is_matrix
            self.data  = np.zeros((rows, cols))
        else:
            self.is_matrix = None
            self.data = None
    
    def init_T(self):
        self.is_matrix = False
        self.data  = np.ones(m+1)
        self.data[0] = 0

    def init_T1(self):
        self.is_matrix = True
        self.data  = np.ones(m+1)
        self.data[0] = 0

    def mul_matrix_by_vector(ps_matrix, ps_vector):
        new_rows = min(ps_matrix.shape[0] + ps_vector.shape[0], MAX_n + 1)
        new_cols = min(ps_matrix.shape[1], MAX_w + 1)

        result = PowerSeries(new_rows, new_cols)

        for matrix_i in range(min(ps_matrix.shape[0], new_rows)):
            for matrix_j in range(min(ps_matrix.shape[1], new_cols)):
                
                matrix_coeff = ps_matrix[matrix_i, matrix_j]

                for other_i in range(min(ps_vector.shape[0], new_rows-matrix_i)):
                    result.data[matrix_i+other_i, matrix_j] += matrix_coeff*ps_vector[other_i]

        return result
    
    def __mul__(self, other):
        ps_self = self.data if (not self.is_matrix or len(self.data.shape) > 1) else np.diag(self.data)
        ps_other = other.data if (not other.is_matrix or len(other.data.shape) > 1) else np.diag(other.data)

        # if (not self.is_matrix or len(self.data.shape) > 1):
        #     ps_self = self.data
        # else: 
        #     ps_self = np.diag(self.data)
        #     ps_self[0,0] = 0
        # if (not other.is_matrix or len(other.data.shape) > 1):
        #     ps_other = other.data
        # else: 
        #     ps_other = np.diag(other.data)
        #     ps_other[0,0] = 0

        if (self.is_matrix and other.is_matrix):
            new_rows = min(ps_self.shape[0] + ps_other.shape[0], MAX_n + 1)
            new_cols = min(ps_self.shape[1] + ps_other.shape[1], MAX_w + 1)

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
            #                 result.data[i,j] += self.data[k,l]*other.data[i-k,j-l]
            
            #option 2:
            for self_i in range(min(ps_self.shape[0], new_rows)):
                for self_j in range(min(ps_self.shape[1], new_cols)):
                    
                    self_coeff = ps_self[self_i, self_j]

                    for other_i in range(min(ps_other.shape[0], new_rows-self_i)):
                        for other_j in range(min(ps_other.shape[1], new_cols-self_j)):
                            result.data[self_i+other_i, self_j+other_j] += self_coeff*ps_other[other_i, other_j]
            
        elif(self.is_matrix):
            return PowerSeries.mul_matrix_by_vector(ps_self, ps_other)

        elif(other.is_matrix):
            return PowerSeries.mul_matrix_by_vector(ps_other, ps_self)
        
        else:  # both vectors
            new_rows = min(ps_self.shape[0] + ps_other.shape[0], MAX_n + 1)
            new_cols = 1

            result = PowerSeries(new_rows, new_cols, is_matrix=False)
            
            for self_i in range(min(ps_self.shape[0], new_rows)):

                self_coeff = ps_self[self_i]

                for other_i in range(min(ps_other.shape[0], new_rows-self_i)):
                    result.data[self_i+other_i] += self_coeff*ps_other[other_i]
    
        return result
    
    def add_vector_to_self(self, ps_vector):

        # note: self.is_matrix must be true
        # note: ps_other shape should be in this point smaller or equal to self

        for other_i in range(ps_vector.shape[0]):
            self.data[other_i,0] += ps_vector[other_i]

        return self
    
    def add_matrix_to_self(self, ps_matrix):

        # note: self.is_matrix must be true
        # note: ps_other shape should be in this point smaller or equal to self

        for matrix_i in range(ps_matrix.shape[0]):
            for matrix_j in range(ps_matrix.shape[1]):
                self.data[matrix_i, matrix_j] += ps_matrix[matrix_i, matrix_j]
        
        return self
    
    def __add__(self, other):
        ps_self = self.data if (not self.is_matrix or len(self.data.shape) > 1) else np.diag(self.data)
        ps_other = other.data if (not other.is_matrix or len(other.data.shape) > 1) else np.diag(other.data)
        
        if (self.is_matrix and other.is_matrix):
            new_rows = max(ps_self.shape[0], ps_other.shape[0])
            new_cols = max(ps_self.shape[1], ps_other.shape[1])

            result = PowerSeries(new_rows, new_cols)

            result.add_matrix_to_self(ps_self)

            return result.add_matrix_to_self(ps_other)
        
        elif (self.is_matrix):
            new_rows = max(ps_self.shape[0], ps_other.shape[0])
            new_cols = ps_self.shape[1]

            result = PowerSeries(new_rows, new_cols)

            result.add_matrix_to_self(ps_self)

            return result.add_vector_to_self(ps_other)

        elif (other.is_matrix):
            new_rows = max(ps_self.shape[0], ps_other.shape[0])
            new_cols = ps_other.shape[1]

            result = PowerSeries(new_rows, new_cols)

            result.add_vector_to_self(ps_self)

            return result.add_matrix_to_self(ps_other)

        # else - both vectors
        new_rows = max(ps_self.shape[0], ps_other.shape[0])

        result = PowerSeries(new_rows, 1, False)

        result.add_vector_to_self(ps_self)

        return result.add_vector_to_self(ps_other)

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