# GPD log likelihood

GPD log likelihood

## Usage

``` r
GPD_log_likelihood(
  Y,
  threshold,
  scale,
  shape,
  obs_weights = NULL,
  ill_defined_value = -10^6
)
```

## Arguments

- Y:

  Data observations.

- threshold:

  GPD threshold value.

- scale:

  Scale parameter.

- shape:

  Shape parameter.

- obs_weights:

  Optional observation weights for weighted likelihood.

- ill_defined_value:

  Value to return if the arguments are out of support (e.g. negative
  scale, or non-positive arguments to logarithms).

## Value

The GPD log-likelihood evaluated at the given parameters, given the data
observations `Y`.
