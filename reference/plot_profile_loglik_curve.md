# Profile likelihood curve plot

Plot the profile likelihood curve for the desired parameter with the
confidence line.

## Usage

``` r
plot_profile_loglik_curve(
  profile_fct_object,
  prop_below = NULL,
  legend.position = "bottom"
)
```

## Arguments

- profile_fct_object:

  Object returned by either
  [`GPD_profile_loglik_curve()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_loglik_curve.md)
  or
  [`GEV_profile_loglik_curve()`](https://opasche.github.io/ExtremeCI/reference/GEV_profile_loglik_curve.md).

- prop_below:

  (Optional) Proportion of the distance below the confidence line to
  show in the plot, to crop y-axis.

- legend.position:

  Position of the legend to the side of the plot. Can be one of
  `"bottom"` (default), `"right"`, `"top"`, `"left"`, or `"none"`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object showing the profile likelihood curve for the desired parameter
with the confidence line.
