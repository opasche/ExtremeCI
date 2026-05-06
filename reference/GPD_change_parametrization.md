# GPD parameter vector reparametrization

GPD parameter vector reparametrization

## Usage

``` r
GPD_change_parametrization(
  parameters,
  threshold = 0,
  threshold_lvl = 0,
  parametrization = c("classical", "orthogonal", "quantile", "endpoint"),
  new_parametrization = c("classical", "orthogonal", "quantile", "endpoint"),
  quantile_lvl = 1 - (1/100),
  nbsca = 1,
  nbsha = 1
)
```

## Arguments

- parameters:

  Vector of GPD parameters in the internal format.

- threshold:

  GPD threshold value.

- threshold_lvl:

  Probability level of the threshold `threshold`.

- parametrization:

  Current parametrization of `parameters`.

- new_parametrization:

  Desired new parametrization.

- quantile_lvl:

  Quantile probability level for the `'quantile'` parameter.

- nbsca:

  Number of scale parameter coefficients (i.e. one plus the number of
  scale parameter covariates).

- nbsha:

  Number of shape parameter coefficients (i.e. one plus the number of
  shape parameter covariates).

## Value

The vector of GPD parameters reparametrized from `parametrization` to
`new_parametrization`.
