# Internal function for the GEV profile likelihood

Internal function for the GEV profile likelihood

## Usage

``` r
GEV_profile_loglik_internal(
  val,
  Z,
  parameter = c("shape", "location", "scale", "return_level", "endpoint"),
  return_period = 100,
  locmat = as.matrix(1),
  scamat = as.matrix(1),
  shamat = as.matrix(1),
  subparam_id = 0,
  orthogonal = FALSE,
  init = NULL,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  ...
)
```

## Arguments

- val:

  .

- Z:

  .

- parameter:

  .

- return_period:

  .

- locmat:

  .

- scamat:

  .

- shamat:

  .

- subparam_id:

  .

- orthogonal:

  .

- init:

  .

- hessian:

  .

- maxit:

  .

- method:

  .

- ...:

  .
