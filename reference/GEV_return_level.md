# Compute return level from GEV parameters

Compute return level from GEV parameters

## Usage

``` r
GEV_return_level(loc = 0, scale = 1, shape, return_period)
```

## Arguments

- loc:

  Location parameter.

- scale:

  Scale parameter.

- shape:

  Shape parameter.

- return_period:

  Return period for the desired return level.

## Value

The return level of the specified GEV distribution with a return period
of `return_period`. In other terms, the quantile of the GEV distribution
at probability level `1 - 1/return_period`.

## Examples

``` r
GEV_return_level(loc=0, scale=1, shape=0.1, return_period=100)
#> [1] 5.840976
```
