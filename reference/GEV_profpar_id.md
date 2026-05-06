# GEV profile parameter vector id

GEV profile parameter vector id

## Usage

``` r
GEV_profpar_id(
  parameter = c("shape", "location", "scale", "return_level", "endpoint"),
  nbloc = 1,
  nbsca = 1,
  nbsha = 1,
  subparam_id = 0
)
```

## Arguments

- parameter:

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

- subparam_id:

  .

## Value

The id of the desired GEV parameter in the specified internal vector
format.
