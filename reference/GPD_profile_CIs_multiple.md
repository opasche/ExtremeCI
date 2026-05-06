# Multi-value conditional GPD profile likelihood confidence intervals

For non-stationary models, the quantile reparametrization depends on
covariate values. This function repeats the profile likelihood procedure
for several covariate values. It enables obtaining a quantile (or
endpoint) curve, with profile-likelihood confidence bands, as a function
of the covariate values.

## Usage

``` r
GPD_profile_CIs_multiple(
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parameter = c("quantile", "endpoint"),
  alpha = 0.05,
  quantile_lvl = 1 - (1/100),
  orthogonal = FALSE,
  X = NULL,
  X_rlvl = NULL,
  scale_cols = NULL,
  shape_cols = NULL,
  init_step_pos = 100,
  init_step_neg = 10,
  tol = 0.01,
  steps_beyond_conf = 5,
  initial_MLE_para = c("classical", "same"),
  max_steps = 10000,
  obs_weights = NULL,
  ill_defined_value = -10^6,
  hessian = TRUE,
  maxit = 1e+06,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  method_prof = c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  parallel_strat = c("none", "multisession", "sequential", "multicore"),
  n_workers = NULL,
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

  Parameter for which to compute the profile likelihood confidence
  intervals.

- alpha:

  Confidence alpha for the profile likelihood confidence intervals.

- quantile_lvl:

  Quantile probability level for the `'quantile'` parameter.

- orthogonal:

  DEPRECATED.

- X:

  Covariate matrix (for conditional/non-stationary fits). Columns should
  be variables, and rows should be observations matching `Y`.

- X_rlvl:

  Covariate matrix at which to reparametrize for the `'quantile'` or
  `'endpoint'` parametrizations (for conditional/non-stationary fits).
  Columns should be variables, and each row should give one covariate
  realization at which to reparametrize and obtain a CI.

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

The GPD profile log-likelihood `(1-alpha)` confidence intervals for the
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

- quantile_lvl:

  Quantile probability level for `'quantile'` parameter (only if
  `parameter==quantile`).
