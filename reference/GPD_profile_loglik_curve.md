# GPD profile log-likelihood curve

GPD profile log-likelihood curve

## Usage

``` r
GPD_profile_loglik_curve(
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parameter = c("shape", "scale", "quantile", "endpoint"),
  subparam_id = 0,
  alpha = 0.05,
  quantile_lvl = 1 - (1/100),
  orthogonal = FALSE,
  X = NULL,
  x_rlvl = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  warmstart_table = NULL,
  stepsize = 0.1,
  steps_beyond_conf = 5,
  initial_MLE_para = c("classical", "same"),
  max_steps = 10000,
  obs_weights = NULL,
  ill_defined_value = -10^6,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  method_prof = c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  ...
)
```

## Arguments

- Y:

  Data observations.

- threshold:

  GPD threshold value.

- threshold_lvl:

  Probability level of the threshold `threshold`.

- parameter:

  Parameter for which to compute the profile log-likelihood.

- subparam_id:

  Index of the parameter coefficient for which to compute the profile
  log-likelihood (for conditional/non-stationary fits).

- alpha:

  Confidence alpha for the profile log-likelihood confidence interval
  (i.e. for the confidence line on the profile plot).

- quantile_lvl:

  Quantile probability level for the `'quantile'` parameter.

- orthogonal:

  DEPRECATED.

- X:

  Covariate matrix (for conditional/non-stationary fits). Columns should
  be variables, and rows should be observations matching `Y`.

- x_rlvl:

  Covariate vector at which to reparametrize for the `'quantile'` or
  `'endpoint'` parametrizations (for conditional/non-stationary fits).

- scale_cols:

  Column indices of `X` to use as covariate for the (conditional) scale
  parameter (for conditional/non-stationary fits).

- shape_cols:

  Column indices of `X` to use as covariate for the (conditional) shape
  parameter (for conditional/non-stationary fits).

- warmstart_table:

  Evaluation table from a previous run.

- stepsize:

  Numerical size of each evaluation step, in the profile parameter's
  scale.

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

- obs_weights:

  Optional observation weights for weighted likelihood.

- ill_defined_value:

  Value to return if the arguments are out of support (e.g. negative
  scale, or non-positive arguments to logarithms).

- hessian:

  Logical. Should a numerically differentiated Hessian matrix be
  returned? See [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  for more details.

- maxit:

  The maximum number of iterations. See
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for more
  details.

- method:

  The optimisation method to be used for the initial maximum likelihood
  optimisation. See
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for more
  details.

- method_prof:

  The optimisation method to be used for the profile likelihood
  optimisation. See
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for more
  details.

- ...:

  Other arguments passed to the `control` argument of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

The GPD profile log-likelihood curve for the desired `parameter`, with
confidence line and resulting `(1-alpha)` confidence interval, as a
`GPD_profileLogLik` object containing:

- mle:

  The estimated maximum likelihood GPD parameters, as a named vector
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

  Index of the profile parameter in the GPD parameter vector.

## References

Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
Values*. Springer.
[doi:10.1007/978-1-4471-3675-0](https://doi.org/10.1007/978-1-4471-3675-0)
.
