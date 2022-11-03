# AlmostBalancedCodeWithoutHomopolymers

This repository is the accompanying code for the project "Almost-Balanced and Maximum Homopolymer-Run Restricted Codes for Data Storage in DNA". For more information about the problem and proposed solutions, please refer to the PDF file of the project report.


## Contents

* **BaselineComputations**: a directoy for computations of optimal bounds.
    * **RLL**: calculates optimal results for the RLL constraint.
        * `RLL3BaselineRedundancy.py`
    * **AlmostBalanced**: calculates optimal results for the imbalance constraint.
        * `AlmostBalancedBaselineRedundancy.py`
    * **AlmostBalancedAndRLL**: calculates optimal results for combined constraints.
        * `PowerSeries.py` and `N4_computation.py`: our custom implementation of the *N4(m, w, n)* formula [by Immink and Cai](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9032172).
        * `calculate_with_backtracking.py`: calculates optimal result with backtracking. Used to validate results of the *N4(m, w, n)* formula implementation for small lengths.
    * **Results**: contains CSV files of the results.
    * `baseline_computations.py`: a script for running all optimal bounds calculations and saving them in seperate CSV files.
    * `strand_requirements.py`: common variables for running scripts and computations.

* **EncodingAlgorithms**: a directoy for computations of results for  and baseline solution and proposed improved solutions.
    * **Results**: contains CSV files of the results.
    * `RLL3andKnuth.py`: a calculation of results for baseline solution.
    * `RLL3andKnuthwithSkips.py`: a calculation of results for the first solution.
    * `RLL3andKnuthbyBlocks.py`: a calculation of results for the second solution.
    * `encoding_computations.py`: a script for running above three calculations and saving them in seperate CSV files.
    * `RLL6andencodings.py`: a calculation of results for the third solution.
    * `common.py`: common variables for running scripts and computations.


## Running the scripts

The project is written in Python 3.

For installing the project's requirements, run from root folder:
```
pip install -r requirements.txt
```
Choose the needed script or implementation, and modify the `main` lines if needed before running the script.


## Authors
* Gili Doweck
* Victoria Goldin
