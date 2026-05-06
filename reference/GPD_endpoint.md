# GPD endpoint

GPD endpoint

## Usage

``` r
GPD_endpoint(
  threshold = 0,
  scale = 1,
  shape,
  ...,
  endpoint_type = c("unrestricted", "upper")
)
```

## Arguments

- threshold:

  GPD threshold value.

- scale:

  Scale parameter.

- shape:

  Shape parameter.

- ...:

  Other optional parameters for the internal endpoint function (Unused).

- endpoint_type:

  Whether to compute the upper endpoint (returns `Inf` if `shape>=0`),
  or the unrestricted endpoint (returns the lower GPD endpoint if
  `shape>=0`).

## Value

The endpoint of the specified GPD distribution.

## Examples

``` r
GPD_endpoint(threshold=0, scale=1, shape=-0.1, endpoint_type='upper')
#> [1] 10
```
