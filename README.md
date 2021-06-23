# Phospho-MV-simulations
Evaluation of differential abundance analysis strategies for phosphoproteomics data with missing values.

We evaluated different testing strategies (i.e. ignoring missing values, imputation, compound hypothesis tests) under 3,234 data scenarios of various levels of sample size, fold change, missing ratio and MNAR ratio. Each scenario was simulated for 500 rounds and a total of 16 methods were compared.

## Requirements
### R (3.6+)
* invgamma
* Bioconductor
* limma
* SDAMS
* pcaMethods
* trust
* ggplot2
* gridExtra
* data.table

### Python (3.7+)
* NumPy (1.19.2+)
* Pandas (1.2.1+)
* Scikit-Learn (0.23.2+)
* SciPy (1.5.2+)

## Simulation
The codes of simulation study are in simulation_HPC.R

Note:
1. The method of SDA was modified by cancelling the data-cleaning step of removing features with low numbers of non-zeros.
2. bPCA imputation was modified by lower the _tol_ parameter in _solve_ function to prevent the error of computationally singular. 

## Results & Plotting
The codes of summarizing adjusted _p_-values of different methods are in results_HPC.py

The codes of plotting figures are in results_plot.R

Only the results data for generating the figures is in this repository due to the limit of file size.
The results data of all methods under all scenarios in the 500-round simulation is available upon request. 
If you have problems, please contact jing.li@sjtu.edu.cn
