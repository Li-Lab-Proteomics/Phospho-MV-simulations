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
The codes of simulation study are in simulation.R

Note:
1. The method of SDA was modified by cancelling the data-cleaning step of removing features with low numbers of non-zeros.
2. bPCA imputation was modified by lower the _tol_ parameter in _solve_ function to prevent from the error of computationally singular. 

## Results & Plotting


