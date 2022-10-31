from math import ceil, floor

REDUNDANCY_LIMIT_PERCENTAGE = 8
REDUNDANCY_LIMIT_LETTERS = 0
IMBALANCE_PERCENTAGE = 5 # The acceptable imbalance of each strand in percentages
m = 3 # Longest homopolymer run
MAX_n_binary = 1000 # Longest binary string length
MAX_n_quaternary = ceil((MAX_n_binary/2.0) + ((0.01*REDUNDANCY_LIMIT_PERCENTAGE/2)*MAX_n_binary)) +\
                   REDUNDANCY_LIMIT_LETTERS # Longest strand length
MAX_w = floor((0.5 + (0.01*IMBALANCE_PERCENTAGE))*MAX_n_quaternary) # Maximum weight of the longest strand
                                                         # (used to restrict unnecessary calculations)

def min_weight(n):
    """ Returns the minimum weight possible for a strand of length n,
        according to the imbalanced percentage defined in this file """

    pure_imbalance = (0.5 - (0.01*IMBALANCE_PERCENTAGE))
    return ceil(pure_imbalance*n)

def max_weight(n):
    """ Returns the maximum weight possible for a strand of length n,
        according to the imbalanced percentage defined in this file """

    pure_imbalance = (0.5 + (0.01*IMBALANCE_PERCENTAGE))
    return floor(pure_imbalance*n)

def min_max_weight(n):
    """ Returns both min and max weight possible for a strand of length n,
        according to calculation in min_weight and max_weight functions.
        This functions also deals with the special case where max_w<min_w """

    min_w = min_weight(n)
    max_w = max_weight(n)
    max_w = max_w + 1 if (max_w<min_w) else max_w

    return min_w, max_w