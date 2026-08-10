# hawkesLib

Package containing various functions for working with Hawkes processes 
(univariate + marked), developed during my STAT489 honours thesis under the 
supervision of Prof. Richard Arnold in the School of Mathematics and Statistics 
at Victoria University of Wellington. 

The main script (`marked-hawkes-npb.R`) contains an implementation of a 
Metropolis-within-Gibbs MCMC sampler to fit a non-parametric kernel to marked 
event time data. Dirichlet process priors are used to model non-parametrically, 
with a step function or a piecewise-linear function implemented. In both cases, 
it is assumed that the kernel is non-increasing (typical for Hawkes models).

There are further scripts, detailed below:
- `sim-functions.R` contains simulation routines for nonhomogeneous Poisson 
processes, univariate Hawkes processes, and marked Hawkes processes. 
- `sim-helpers.R`contains wrappers around the simulation routines and the MCMC 
sampler functions to perform a simulation study, computing bias, RMSE, credible 
interval coverage, and plotting 'true' kernels against posterior estimates
- `plot-utils.R` contains some functions to generate commonly used plots
- `exp-hp-mle.R` contains functions to perform MLE for univariate Hawkes process
models with a parametric exponential kernel.