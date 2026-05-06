# GEV parameter vector reparametrization

GEV parameter vector reparametrization

## Usage

``` r
GEV_change_parametrization(
  parameters,
  parametrization = c("classical", "return_level", "endpoint"),
  new_parametrization = c("classical", "return_level", "endpoint"),
  return_period = 100,
  nbloc = 1,
  nbsca = 1,
  nbsha = 1
)
```

## Arguments

- parameters:

  Vector of GEV parameters in the internal format.

- parametrization:

  Current parametrization of `parameters`.

- new_parametrization:

  Desired new parametrization.

- return_period:

  Return period for the `'return_level'` parameter.

- nbloc:

  Number of location parameter coefficients (i.e. one plus the number of
  shape parameter covariates).

- nbsca:

  Number of scale parameter coefficients (i.e. one plus the number of
  scale parameter covariates).

- nbsha:

  Number of shape parameter coefficients (i.e. one plus the number of
  shape parameter covariates).

## Value

The vector of GEV parameters reparametrized from `parametrization` to
`new_parametrization`.
