---
editor_options: 
  markdown: 
    wrap: 72
---

## `robcdcc` a package for robust cDCC estimation

The `robcdcc` is a `R` package for conditional variance and covariance
estimation, using the robust corrected dynamic conditional correlation
(cDCC) model proposed by [Boudt, Daníelsson and Laurent (2013)](https://doi.org/10.1016/j.ijforecast.2012.06.003) , named BIP-cDCC. The
estimation procedures for the traditional cDCC model are also available.

# Instalation

The development version of `robcdcc` can be installed from
[GitHub](https://github.com/) using:

    library(devtools)
    install_github("EduardoGPinheiro/robcdcc")

# Functions

Some useful functions in this library:

-   `simCDCC`: Simulate GARCH(1,1)-cDCC(1,1) bivariate returns;

-   `estimateCDCC`: Estimate GARCH(1,1)-cDCC(1,1) parameters;

-   `robust_estimateCDCC`: Estimate BIP-GARCH(1,1)-BIP-cDCC(1,1)
    parameters;
    
-   `calc_ht`: Calculate estimated conditional variance for GARCH(1,1);

-   `robust_calc_ht`: Calculate estimated conditional variance for
    BIP-GARCH(1,1);
    
-   `calc_Rt`: Calculate estimated conditional correlation matrix for
    cDCC(1,1);
    
-   `robust_calc_Rt`: Calculate robust estimated conditional correlation matrix
    for BIP-cDCC(1,1).

# Examples

We'll illustrate `robcdcc` by a simulated example. Let $\mathbf{r}_t$ be
a bivariate series of returns generated by a GARCH(1,1)-cDCC(1,1)
process given by:

$$
\mathbf{r}_t = \mathbf{H}_t^{1/2} \boldsymbol{\epsilon}_t, \quad \text{given} 
\quad \boldsymbol{\epsilon}_t \sim \mathcal{N}(0, I_2)
$$

where

$$ 
h_{t, 1} = 0,10 + 0,10 r_{t-1, 1}^2 + 0,80 h_{t-1, 1}, \\ 
h_{t, 2} = 0,10 + 0,20 r_{t-1, 2}^2 + 0,70 h_{t-1, 2}, \\
\mathbf{Q}_t = (1 - 0.1 - 0.8) \mathbf{S} + 0.1 \mathbf{Q}_{t-1}^{1/2
*} \mathbf{r}_{t-1} \mathbf{r}_{t-1}^\top \mathbf{Q}_{t-1}^{1/2
*} + 0.8 \mathbf{Q}_t,\\
\mathbf{R}_t = \mathbf{Q}_{t-1}^{-1/2*}\mathbf{Q}_t\mathbf{Q}_{t-1}^{-1/2}, \\
\mathbf{H}_t = \mathbf{D}_t^{1/2} \mathbf{R}_t \mathbf{D}_t^{1/2}
$$

for $\mathbf{D}_t = diag(h_{t, 1}, h_{t, 2})$ and $\mathbf{S}_{1,2} = 0,4$.

```
# Simulation
S = diag(2) * .4
eta = matrix(c(.1, .1, .8, .1, .2, .7), byrow = TRUE, ncol = 3)
phi = c(.1, .8)

sim_lst = simCDCC(phi=phi, eta=eta, S=S, nobs=1000, ndim=2, seed=1)
rt = sim_lst$rt # returns matrix
ht = sim_lst$ht # conditional variance matrix 
Rt = sim_lst$Rt # conditional correlations matrices

# Estimation
q_results = estimateCDCC(rt=rt) # traditional estimation procedure
r_results = robust_estimateCDCC(rt=rt) # robust estamation procedure
```

