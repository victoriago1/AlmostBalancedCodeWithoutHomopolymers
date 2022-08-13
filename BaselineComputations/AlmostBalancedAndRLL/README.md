
# Baseline computation for both Almost balanced and RLL
This folder contains everything needed to compute the number of strands that follow the requirements as defined in
`strand_requirements.py` (More details on this file are bellow).

Must of the files are relevant to the computation of the function N_4 as it is defined in the article:\n
"Properties and Constructions of Constrained Codes for DNA-Based Data Storage"[^1].
The function N_4, defined in this article, calculates the number of strand of length n and wight w using recursive 
functions and power series.
Please read the article for more information.

This folder contains or uses the following files:
- `../strand_requirements.py`:
    This file holds globally defined parameters representing the requirements from a strand, including -
    - percentage of imbalanced allowed
    - longest length of homopolymer run allowed
    - restrictions for the computations, such as the maximum length of a strand the computation should include

- `calculate_with_backtracking.py`:
    This file computes the number of strands that follow the requirements by brute force.
    This is in order to verify that the computations are correct.

- `PowerSeries.py`:
    This files implements a class "PowerSeries" that is used to easily define the power series defined in the article
    and easily compute the multiplication and addition of different power series.

- `N4_computation.py`:
    This file computes the function N_4 defined in the article and saves the results under the file '../Results'.
    Make sure to run this from the git root folder

## requirements:
#TODO - fix?
```bash
pip install math
pip install numpy
pip install pandas
pip install tqdm
```

[^1]: K. A. S. Immink and K. Cai, "Properties and Constructions of Constrained Codes for DNA-Based Data Storage," in IEEE Access, vol. 8, pp. 49523-49531, 2020, doi: 10.1109/ACCESS.2020.2980036.