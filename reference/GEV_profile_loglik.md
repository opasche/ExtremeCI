# GEV profile log-likelihood

GEV profile log-likelihood

## Usage

``` r
GEV_profile_loglik(
  val,
  Z,
  parameter = c("shape", "location", "scale", "return_level", "endpoint"),
  subparam_id = 0,
  return_period = 100,
  orthogonal = FALSE,
  X = NULL,
  x_rlvl = NULL,
  loc_cols = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  init = NULL,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  ...
)
```

## Arguments

- val:

  Parameter value at which to evaluate the GEV profile log-likelihood.

- Z:

  Block maxima observations.

- parameter:

  Parameter for which to compute the profile log-likelihood.

- subparam_id:

  Index of the parameter coefficient for which to compute the profile
  log-likelihood (for conditional/non-stationary fits).

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

The GEV profile log-likelihood of `parameter` evaluated at `val`, given
the data, as a `GEV_profML` object containing:

- param_val:

  (Named) parameter value at which the GEV profile log-likelihood was
  evaluated.

- param_name:

  Name of the evaluated profile likelihood parameter.

- mle_other:

  Maximum-likelihood estimate of the other GEV parameters.

- loglik:

  Profile GEV log-likelihood value of `parameter` evaluated at `val`,
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
