# Internal function for the GPD profile likelihood

Internal function for the GPD profile likelihood

## Usage

``` r
GPD_profile_loglik_internal(
  val,
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parameter = c("shape", "scale", "quantile", "endpoint"),
  quantile_lvl = 1 - (1/100),
  scamat = as.matrix(1),
  shamat = as.matrix(1),
  subparam_id = 0,
  orthogonal = FALSE,
  obs_weights = NULL,
  ill_defined_value = -10^6,
  init = NULL,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  verbose = 1,
  ...
)
```

## Arguments

- val:

  .

- Y:

  .

- threshold:

  .

- threshold_lvl:

  .

- parameter:

  .

- quantile_lvl:

  .

- scamat:

  .

- shamat:

  .

- subparam_id:

  .

- orthogonal:

  .

- obs_weights:

  Optional observation weights for weighted likelihood.

- ill_defined_value:

  Value to return if the arguments are out of support (e.g. negative
  scale, or non-positive arguments to logarithms).

- init:

  .

- hessian:

  .

- maxit:

  .

- method:

  .

- verbose:

  Verbose level, as integer.

- ...:

  .
