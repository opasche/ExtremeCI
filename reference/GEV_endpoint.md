# GEV endpoint

GEV endpoint

## Usage

``` r
GEV_endpoint(
  loc = 0,
  scale = 1,
  shape,
  ...,
  endpoint_type = c("unrestricted", "upper", "lower")
)
```

## Arguments

- loc:

  Location parameter.

- scale:

  Scale parameter.

- shape:

  Shape parameter.

- ...:

  Other optional parameters for the internal endpoint function (Unused).

- endpoint_type:

  Whether to compute the upper endpoint (returns `Inf` if `shape>=0`),
  or the unrestricted endpoint (returns the lower GEV endpoint if
  `shape>=0`).

## Value

The endpoint of the specified GEV distribution.

## Examples

``` r
GEV_endpoint(loc=0, scale=1, shape=-0.1, endpoint_type='upper')
#> [1] 10
```
