from math import ceil, floor

"""
Calculates by backtracking the number of strands that hold the homopolymer constraint
and the GC/AT 5% almost balanced constraint.
"""

SIGMA = ['a','t','c','g']


def calc_longest_homopolymer(str) -> int:

    """ Returns the maximum homopolymer length. """

    if not str:
        return 0
    
    longest = 1
    prev = None
    current_len = 0
    for char in str:
        if char == prev:
            current_len += 1
            if current_len > longest:
                longest = current_len
        else:
            prev = char
            current_len = 1
    return longest


def calc_weight(str):

    """ Returns the wheight of the strand. """

    return str.count('a') + str.count('t')


def calc_strands(n, m, verbose=False) -> int:

    """ Returns the number of strands of length n and maximum homopolymer length of m.
        If verbose==True, all of the strands that hold the conditions will be printed as well. """
    
    w_min = ceil(0.45 * n)
    w_max = floor(0.55 * n)
    if w_max < w_min:
        w_min = w_max = n//2  # by definition of balanced in Knuth
    return _calc_strands_wrapper(n=n, m=m, w_min=w_min, w_max=w_max, str='', verbose=verbose)


def _calc_strands_wrapper(n, m, w_min, w_max, str, verbose) -> int:
    if calc_longest_homopolymer(str) > m:
        return 0

    if len(str) == n:
        weight = calc_weight(str)
        if weight >= w_min and weight <= w_max:
            if verbose:
                print(str)
            return 1
        return 0
    
    counter = 0
    for char in SIGMA:
        counter += _calc_strands_wrapper(n=n, m=m, w_min=w_min, w_max=w_max, str=str + char, verbose=verbose)
    return counter


if __name__ == "__main__":
    # print(calc_strands(n=4, m=3, verbose=True))

    m = 3
    for n in range(2, 13):
        print('Number of strands for n={} and m={}: {}'.format(n, m, calc_strands(n=n, m=m, verbose=False)))
    