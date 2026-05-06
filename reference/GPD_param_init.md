# Initial GPD parameter vector defaults for profile optimization

Initial GPD parameter vector defaults for profile optimization

## Usage

``` r
GPD_param_init(
  Y,
  threshold = 0,
  threshold_lvl = 0,
  nbsca,
  nbsha,
  parametrization = c("classical", "orthogonal", "quantile", "endpoint"),
  quantile_lvl,
  obs_weights = NULL
)
```

## Arguments

- Y:

  .

- threshold:

  .

- threshold_lvl:

  .

- nbsca:

  .

- nbsha:

  .

- parametrization:

  .

- quantile_lvl:

  .

## Value

The initial parameter values as a vector, in the correct internal
format.
