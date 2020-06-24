# cenzus
Multiple Systems Estimation

The package is developed for the analysis of multiple population registers each containing the self-reported covariate ethinicity, which includes missing values and potentially measurement error. The function `mse` fits a loglinear model with  optional latent classes to account for measurement error, and returns the population size estimates. The function `boot_mse` performs a semi-parametric bootstrap for the computation of confidence intervals.
