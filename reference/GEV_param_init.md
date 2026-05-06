# Initial GEV parameter vector defaults for profile optimization

Initial GEV parameter vector defaults for profile optimization

## Usage

``` r
GEV_param_init(
  Z,
  nbloc,
  nbsca,
  nbsha,
  parametrization = c("classical", "return_level", "endpoint"),
  return_period
)
```

## Arguments

- Z:

  .

- nbloc:

  Number of location parameter coefficients (i.e. one plus the number of
  shape parameter covariates).

- nbsca:

  Number of scale parameter coefficients (i.e. one plus the number of
  scale parameter covariates).

- nbsha:

  Number of shape parameter coefficients (i.e. one plus the number of
  shape parameter covariates).

- parametrization:

  .

- return_period:

  .

## Value

The initial parameter values as a vector, in the correct internal
format.
