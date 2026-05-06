# Multi-value conditional GEV profile likelihood confidence intervals

For non-stationary models, the return level reparametrization depends on
covariate values. This function repeats the profile likelihood procedure
for several covariate values. It enables obtaining a return-level (or
endpoint) curve, with profile-likelihood confidence bands, as a function
of the covariate values.

## Usage

``` r
GEV_profile_CIs_multiple(
  Z,
  parameter = c("return_level", "endpoint"),
  alpha = 0.05,
  return_period = 100,
  orthogonal = FALSE,
  X = NULL,
  X_rlvl = NULL,
  loc_cols = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  init_step_pos = 100,
  init_step_neg = 10,
  tol = 0.01,
  steps_beyond_conf = 5,
  initial_MLE_para = c("classical", "same"),
  max_steps = 10000,
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

- parameter:

  Parameter for which to compute the profile likelihood confidence
  intervals.

- alpha:

  Confidence alpha for the profile likelihood confidence intervals.

- return_period:

  Return period for the `'return_level'` parameter.

- orthogonal:

  DEPRECATED.

- X:

  Covariate matrix (for conditional/non-stationary fits). Columns should
  be variables, and rows should be observations matching `Y`.

- X_rlvl:

  Covariate matrix at which to reparametrize for the `'return_level'` or
  `'endpoint'` parametrizations (for conditional/non-stationary fits).
  Columns should be variables, and each row should give one covariate
  realization at which to reparametrize and obtain a CI.

- loc_cols:

  Column indices of `X` to use as covariate for the (conditional)
  location parameter (for conditional/non-stationary fits).

- scale_cols:

  Column indices of `X` to use as covariate for the (conditional) scale
  parameter (for conditional/non-stationary fits).

- shape_cols:

  Column indices of `X` to use as covariate for the (conditional) shape
  parameter (for conditional/non-stationary fits).

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

- parallel_strat:

  Parallel strategy. One of `"sequential"` (default), `"multisession"`,
  `"multicore"`, or `"mixed"`.

- n_workers:

  A positive numeric scalar or a function specifying the maximum number
  of parallel futures that can be active at the same time before
  blocking. If a function, it is called without arguments when the
  future is created and its value is used to configure the workers. The
  function should return a numeric scalar. Defaults to
  [`future::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)`-1`
  if `NULL` (default), with `"multicore"` constraint in the relevant
  case. Ignored if `strategy=="sequential"`.

- ...:

  Other arguments passed to the `control` argument of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

The GEV profile log-likelihood `(1-alpha)` confidence intervals for the
desired `parameter`, for each desired covariate values, as a
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html),
with columns:

- obs:

  Index of the observation (i.e. row) of X_rlvl for which the parameter
  estimate and CI was computed.

- `<parameter name>`:

  Conditional estimate of the `parameter`.

- ci_down:

  Lower endpoint of the conditional `(1-alpha)` profile-likelihood
  confidence interval for `parameter`.

- ci_up:

  Upper endpoint of the conditional `(1-alpha)` profile-likelihood
  confidence interval for `parameter`.

- parameter:

  Name of the parameter for which the estimates and CIs were computed.

- alpha:

  Confidence alpha for the profile likelihood confidence intervals.

- return_period:

  Return period for the `'return_level'` parameter (only if
  `parameter==return_level`).
