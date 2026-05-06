# Optimisation evaluation step for the GPD profile likelihood

Optimisation evaluation step for the GPD profile likelihood

## Usage

``` r
optim_step_GPD_profile(
  a,
  val,
  Y,
  threshold = 0,
  threshold_lvl = 0,
  parametrization = c("classical", "orthogonal", "quantile", "endpoint"),
  id_param,
  quantile_lvl = 1 - (1/100),
  scamat = as.matrix(1),
  shamat = as.matrix(1),
  negative = TRUE,
  orthogonal = FALSE,
  obs_weights = NULL,
  ill_defined_value = -10^6
)
```

## Arguments

- a:

  .

- val:

  .

- Y:

  .

- threshold:

  .

- threshold_lvl:

  .

- parametrization:

  .

- id_param:

  .

- quantile_lvl:

  .

- scamat:

  .

- shamat:

  .

- negative:

  .

- orthogonal:

  .

- obs_weights:

  Optional observation weights for weighted likelihood.

- ill_defined_value:

  Value to return if the arguments are out of support (e.g. negative
  scale, or non-positive arguments to logarithms).
