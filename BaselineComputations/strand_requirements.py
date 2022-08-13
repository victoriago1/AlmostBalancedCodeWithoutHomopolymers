from math import ceil, floor

m = 3 # Longest homopolymer run
MAX_n = 50 # Longest string length
IMBALANCE_PERCENTAGE = 5 # The acceptable imbalance of each strand in percentages
MAX_w = floor((0.5 + (0.01*IMBALANCE_PERCENTAGE))*MAX_n) # Maximum wight of the longest strand
                                                         # (used to restrict unnecessary calculations)

def min_wight(n):
    """ Returns the minimum wight possible for a strand of length n,
        according to the imbalanced percentage defined in this file """

    pure_imbalance = (0.5 - (0.01*IMBALANCE_PERCENTAGE))
    return ceil(pure_imbalance*n)

def max_wight(n):
    """ Returns the maximum wight possible for a strand of length n,
        according to the imbalanced percentage defined in this file """

    pure_imbalance = (0.5 + (0.01*IMBALANCE_PERCENTAGE))
    return floor(pure_imbalance*n)

def min_man_w(n):
    """ Returns both min and max wight possible for a strand of length n,
        according to calculation in min_wight and max_wight functions.
        This functions also deals with the special case where max_w<min_w """

    min_w = min_wight(n)
    max_w = max_wight(n)
    max_w = max_w + 1 if (max_w<min_w) else max_w

    return min_w, max_w