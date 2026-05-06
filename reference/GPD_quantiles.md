# Compute extreme quantile from GPD parameters

Compute extreme quantile from GPD parameters

## Usage

``` r
GPD_quantiles(quantile_lvl, threshold_lvl, threshold, scale, shape)
```

## Arguments

- quantile_lvl:

  Probability level of the desired extreme quantile.

- threshold_lvl:

  Probability level of the GPD threshold.

- threshold:

  GPD threshold value.

- scale:

  Value(s) for the GPD scale parameter.

- shape:

  Value(s) for the GPD shape parameter.

## Value

The quantile value at probability level `quantile_lvl`.

## Examples

``` r
GPD_quantiles(quantile_lvl=0.999, threshold_lvl=0.95, threshold=0, scale=1, shape=0.1)
#> [1] 4.787576
```
