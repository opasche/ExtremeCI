# Internal GPD log-likelihood function format for optimisation

Internal GPD log-likelihood function format for optimisation

## Usage

``` r
GPD_log_likelihood_optim(
  a,
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parametrization = c("classical", "orthogonal", "quantile", "endpoint"),
  quantile_lvl = 1 - (1/100),
  scamat = as.matrix(1),
  shamat = as.matrix(1),
  negative = FALSE,
  obs_weights = NULL,
  ill_defined_value = -10^6
)
```

## Arguments

- a:

  .

- Y:

  .

- threshold:

  .

- threshold_lvl:

  .

- parametrization:

  .

- quantile_lvl:

  .

- scamat:

  .

- shamat:

  .

- negative:

  .

- obs_weights:

  Optional observation weights for weighted likelihood.

- ill_defined_value:

  Value to return if the arguments are out of support (e.g. negative
  scale, or non-positive arguments to logarithms).
