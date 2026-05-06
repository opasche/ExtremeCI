# GEV_paraboot_CIs_multiple (IN DEVELOPPEMENT)

GEV_paraboot_CIs_multiple (IN DEVELOPPEMENT)

## Usage

``` r
GEV_paraboot_CIs_multiple(
  Z,
  alpha = 0.05,
  return_period = 100,
  R = 1000,
  endpoint_type = c("upper", "lower", "unrestricted"),
  orthogonal = FALSE,
  X = NULL,
  X_rlvl = NULL,
  loc_cols = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  bootstrap_X = FALSE,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  parallel_strat = c("none", "multisession", "sequential", "multicore"),
  n_workers = NULL,
  ...
)
```

## Arguments

- Z:

  Block maxima observations.

- alpha:

  .

- return_period:

  .

- R:

  .

- endpoint_type:

  .

- orthogonal:

  .

- X:

  .

- X_rlvl:

  .

- loc_cols:

  .

- scale_cols:

  .

- shape_cols:

  .

- bootstrap_X:

  .

- hessian:

  .

- maxit:

  .

- method:

  .

- parallel_strat:

  .

- n_workers:

  .

- ...:

  .

## Value

Parametric bootstrap CIs
