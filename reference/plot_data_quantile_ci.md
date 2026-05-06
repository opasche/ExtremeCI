# Non-stationary confidence bands plot

Plot the estimated conditional quantile, return level, or endpoint, with
confidence bands for non-stationary models.

## Usage

``` r
plot_data_quantile_ci(
  quantiles,
  q_down,
  q_up,
  time_index = seq_along(quantiles),
  x_label = "X",
  y_label = "Y",
  Y = NULL,
  event_index = NULL,
  obs_label = "Fitting obs.",
  event_label = "Event",
  legend.position = "bottom"
)
```

## Arguments

- quantiles:

  Vector of estimated quantiles, return levels, or endpoints.

- q_down:

  Vector of lower confidence band values for `quantiles`.

- q_up:

  Vector of upper confidence band values for `quantiles`.

- time_index:

  Vector of time indices corresponding to the `quantiles`. Defaults to
  `seq_along(quantiles)`.

- x_label:

  Label for the x-axis. Defaults to `"X"`.

- y_label:

  Label for the y-axis. Defaults to `"Y"`.

- Y:

  (Optional) Vector of observations to add to the plot and compare to
  the `quantiles.`

- event_index:

  (Optional) Index of an event (observation) to be highlighted.

- obs_label:

  Label for the fitting observations. Defaults to `"Fitting obs."`.

- event_label:

  (Optional) Label for the highlighted event observation. Defaults to
  `"Event"`.

- legend.position:

  Position of the legend to the side of the plot. Can be one of
  `"bottom"` (default), `"right"`, `"top"`, `"left"`, or `"none"`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object showing the estimated conditional quantile with confidence bands.
