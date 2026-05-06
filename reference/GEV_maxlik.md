# Maximum-likelihood GEV estimate

Maximum-likelihood GEV estimate

## Usage

``` r
GEV_maxlik(
  Z,
  parametrization = c("classical", "return_level", "endpoint"),
  return_period = 100,
  orthogonal = FALSE,
  X = NULL,
  x_rlvl = NULL,
  loc_cols = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  out_param = parametrization,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  ...
)
```

## Arguments

- Z:

  Block maxima observations.

- parametrization:

  Likelihood parametrization. Alternatives to `classical` substitute for
  the location parameter.

- return_period:

  Return period for the `'return_level'` parametrization.

- orthogonal:

  DEPRECATED.

- X:

  Covariate matrix (for conditional/non-stationary fits).

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

- out_param:

  Additional output parametrization (same as `parametrization`, by
  default). If `out_param != parametrization`, the parameters are
  reparametrized from `parametrization` to `out_param` after estimation,
  in a separate output.

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

The fitted maximum-likelihood GEV as a `GEV_ML` object, containing:

- mle:

  The estimated maximum likelihood GEV parameters, as a named vector
  (expressed in `parametrization`).

- loglik:

  The log-likelihood of the estimated parameters, given the data

- conv:

  Whether the optimisation procedure converged. See the `convergence`
  output of [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for
  more details.

- hessian:

  The loglikelihood hessian evaluated at the estimated parameters, given
  the data.

- parametrization:

  Likelihood parametrization.

- out_mle:

  The estimated maximum likelihood GEV parameters, reparametrized in
  `out_param`.

- out_parametrization:

  Additional output parametrization.
