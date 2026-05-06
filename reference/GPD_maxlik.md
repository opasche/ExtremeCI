# Maximum-likelihood GPD estimate

Maximum-likelihood GPD estimate

## Usage

``` r
GPD_maxlik(
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parametrization = c("classical", "orthogonal", "quantile", "endpoint"),
  quantile_lvl = 1 - (1/100),
  orthogonal = FALSE,
  X = NULL,
  x_rlvl = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  out_param = parametrization,
  obs_weights = NULL,
  ill_defined_value = -10^6,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  verbose = 1,
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

- parametrization:

  Likelihood parametrization. Alternatives to `classical` substitute for
  the scale parameter.

- quantile_lvl:

  Quantile probability level for the `'quantile'` parametrization.

- orthogonal:

  DEPRECATED.

- X:

  Covariate matrix (for conditional/non-stationary fits).

- x_rlvl:

  Covariate vector at which to reparametrize for the `'quantile'` or
  `'endpoint'` parametrizations (for conditional/non-stationary fits).

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

  The optimisation method to be used. See
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) for more
  details.

- verbose:

  Verbose level, as integer.

- ...:

  Other arguments passed to the `control` argument of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

The fitted maximum-likelihood GPD as a `GPD_ML` object, containing:

- mle:

  The estimated maximum likelihood GPD parameters, as a named vector
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

  The estimated maximum likelihood GPD parameters, reparametrized in
  `out_param`.

- out_parametrization:

  Additional output parametrization.
