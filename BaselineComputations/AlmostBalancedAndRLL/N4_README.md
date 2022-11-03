
# Baseline computation for both Almost balanced and RLL constraints
This folder contains everything needed to compute the number of strands that follow the requirements as defined in
`strand_requirements.py` (More details on this file are bellow).

Most of the files are relevant to the computation of the function *N4(m, w, n)* as it is defined in the article:
"Properties and Constructions of Constrained Codes for DNA-Based Data Storage"[^1].
The function *N4(m, w, n)*, defined in this article, calculates the number of strand of length n and weight w using recursive 
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
    This files implements a custom class "PowerSeries", including initializing the basic power series that are used in the article. The class easily computes multiplication and addition of different power series.

- `N4_computation.py`:
    This file computes the function *N4(m, w, n)* defined in the article and saves the results under the directory '../Results'.
    The script is meant to be run from the repository's root folder.

## Requirements

For installing the project's requirements, run:
```
pip install -r requirements.txt
```

[^1]: K. A. S. Immink and K. Cai, ["Properties and Constructions of Constrained Codes for DNA-Based Data Storage"](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9032172), in IEEE Access, vol. 8, pp. 49523-49531, 2020, doi: 10.1109/ACCESS.2020.2980036.
