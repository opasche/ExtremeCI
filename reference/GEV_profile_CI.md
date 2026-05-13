# GEV profile CI using binary search

GEV profile CI using binary search

## Usage

``` r
GEV_profile_CI(
  Z,
  parameter = c("shape", "location", "scale", "return_level", "endpoint"),
  subparam_id = 0,
  alpha = 0.05,
  return_period = 100,
  orthogonal = FALSE,
  X = NULL,
  x_rlvl = NULL,
  loc_cols = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  warmstart_table = NULL,
  init_step_pos = 100,
  init_step_neg = 10,
  tol = 0.01,
  steps_beyond_conf = 5,
  initial_MLE_para = c("classical", "same"),
  max_steps = 1000,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  verbose = 1,
  ...
)
```

## Arguments

- Z:

  Block maxima observations.

- parameter:

  Parameter for which to compute the profile log-likelihood.

- subparam_id:

  Index of the parameter coefficient for which to compute the profile
  log-likelihood (for conditional/non-stationary fits).

- alpha:

  Confidence alpha for the profile likelihood confidence interval (i.e.
  for the confidence line on the profile plot).

- return_period:

  Return period for the `'return_level'` parameter.

- orthogonal:

  DEPRECATED.

- X:

  Covariate matrix (for conditional/non-stationary fits). Columns should
  be variables, and rows should be observations matching `Y`.

- x_rlvl:

  Covariate vector at which to reparametrize for the `'return_level'` or
  `'endpoint'` parametrizations (for conditional/non-stationary fits).

- loc_cols:

  Column indices of `X` to use as covariate for the (conditional)
  location parameter (for conditional/non-stationary fits).

- scale_cols:

  Column indices of `X` to use as covariate for the (conditional) scale
  parameter (for conditional/non-stationary fits).

- shape_cols:

  Column indices of `X` to use as covariate for the (conditional) shape
  parameter (for conditional/non-stationary fits).

- warmstart_table:

  Evaluation table from a previous run.

- init_step_pos:

  Initial numerical size of each evaluation step to the right, in the
  profile parameter's scale.

- init_step_neg:

  Initial numerical size of each evaluation step to the left, in the
  profile parameter's scale.

- tol:

  Numerical tolerance for convergence, in the profile parameter's scale.

- steps_beyond_conf:

  Number of additional steps to take (in each direction) after the
  profile log-likelihood values reach below the confidence line.

- initial_MLE_para:

  Parametrization used for the initial maximum likelihood estimate
  (defaults to classical, for better stability).

- max_steps:

  Maximum number of steps taken (in each direction). If the confidence
  line was not reached, the corresponding confidence interval endpoint
  will be infinite.

- hessian:

  Logical. Should a numerically differentiated Hessian matrix be
  returned? See [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  for more details.

- maxit:

  The maximum number of iterations. See
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for more
  details.

- method:

  The optimisation method to be used. See
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for more
  details.

- verbose:

  Verbose level, as integer.

- ...:

  Other arguments passed to the `control` argument of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

The GEV profile log-likelihood confidence interval for the desired
`parameter`, with confidence line and resulting `(1-alpha)` confidence
interval, as a `GEV_profileLogLik` object containing:

- mle:

  The estimated maximum likelihood GEV parameters, as a named vector
  (expressed in the profile parametrization).

- ci:

  Length-two vector containing the lower and upper endpoints of the
  desired profile likelihood confidence interval.

- profile_loglik:

  Named matrix containing the profile loglikelihood value (Column 2) for
  each considered profile parameter value (Column 1).

- conf_line:

  Confidence line for the desired profile likelihood confidence
  interval. See e.g. Coles (2001) for more details.

- eval_table:

  Tibble
  ([`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html))
  containing the history of profile log-likelihood evaluation values,
  and related metadata.

- param_name:

  Name of the profiled parameter (infered, for debugging purposes).

- parameter:

  Name of the profiled parameter (given).

- parametrization:

  Parametrization used for the profile likelihood.

- subparam_id:

  Index of the parameter coefficient for which the profile
  log-likelihood was computed.

- id_param:

  Index of the profile parameter in the GEV parameter vector.

## References

Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
Values*. Springer.
[doi:10.1007/978-1-4471-3675-0](https://doi.org/10.1007/978-1-4471-3675-0)
.
