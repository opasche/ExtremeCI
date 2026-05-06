# GPD profile log-likelihood

GPD profile log-likelihood

## Usage

``` r
GPD_profile_loglik(
  val,
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parameter = c("shape", "scale", "quantile", "endpoint"),
  subparam_id = 0,
  quantile_lvl = 1 - (1/100),
  orthogonal = FALSE,
  X = NULL,
  x_rlvl = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  obs_weights = NULL,
  ill_defined_value = -10^6,
  init = NULL,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  ...
)
```

## Arguments

- val:

  Parameter value at which to evaluate the GPD profile log-likelihood.

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

- obs_weights:

  Optional observation weights for weighted likelihood.

- ill_defined_value:

  Value to return if the arguments are out of support (e.g. negative
  scale, or non-positive arguments to logarithms).

- init:

  Optional initial values for the remaining parameter's optimisation
  process, in the correct internal format.

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

- ...:

  Other arguments passed to the `control` argument of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

The GPD profile log-likelihood of `parameter` evaluated at `val`, given
the data, as a `GPD_profML` object containing:

- param_val:

  (Named) parameter value at which the GPD profile log-likelihood was
  evaluated.

- param_name:

  Name of the evaluated profile likelihood parameter.

- mle_other:

  Maximum-likelihood estimate of the other GPD parameters.

- loglik:

  Profile GPD log-likelihood value of `parameter` evaluated at `val`,
  given the data.

- conv:

  Whether the optimisation procedure converged. See the `convergence`
  output of [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for
  more details.

- hessian:

  The loglikelihood hessian evaluated at the estimated parameters, given
  the data.

- parameter:

  Name of the evaluated profile likelihood parameter given as argument
  (redundent).

- parametrization:

  Likelihood parametrization.

- subparam_id:

  Index of the parameter coefficient for which the profile
  log-likelihood was computed.

- id_param:

  Index of the likelihood profile parameter, in the internal parameter
  vector format.
