# Package index

## Profile CIs using peaks over threshold (generalized Pareto distribution)

In most cases, if the data is not aggregated as block-maxima, the
generalized Pareto distribution is the best choice for extrapolating
quantile estimates beyond the range of observations as it makes better
use of available data, at the cost of needing to select an appropriate
threshold value.

- [`GPD_profile_loglik_curve()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_loglik_curve.md)
  : GPD profile log-likelihood curve
- [`GPD_profile_CI()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_CI.md)
  : GPD profile CI using binary search
- [`GPD_profile_CIs_multiple()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_CIs_multiple.md)
  : Multi-value conditional GPD profile likelihood confidence intervals
- [`GPD_profile_loglik()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_loglik.md)
  : GPD profile log-likelihood
- [`GPD_maxlik()`](https://opasche.github.io/ExtremeCI/reference/GPD_maxlik.md)
  : Maximum-likelihood GPD estimate
- [`GPD_quantiles()`](https://opasche.github.io/ExtremeCI/reference/GPD_quantiles.md)
  : Compute extreme quantile from GPD parameters
- [`GPD_endpoint()`](https://opasche.github.io/ExtremeCI/reference/GPD_endpoint.md)
  : GPD endpoint
- [`GPD_log_likelihood()`](https://opasche.github.io/ExtremeCI/reference/GPD_log_likelihood.md)
  : GPD log likelihood
- [`GPD_change_parametrization()`](https://opasche.github.io/ExtremeCI/reference/GPD_change_parametrization.md)
  : GPD parameter vector reparametrization

## Profile CIs using block maxima (generalized extreme value distribution)

When the data at hand is aggregated as block maxima (e.g. yearly
maxima), the generalized extreme value distribution is the natural
choice to extrapolate quantile estimates beyond the range of
observations.

- [`GEV_profile_loglik_curve()`](https://opasche.github.io/ExtremeCI/reference/GEV_profile_loglik_curve.md)
  : GEV profile log-likelihood curve
- [`GEV_profile_CI()`](https://opasche.github.io/ExtremeCI/reference/GEV_profile_CI.md)
  : GEV profile CI using binary search
- [`GEV_profile_CIs_multiple()`](https://opasche.github.io/ExtremeCI/reference/GEV_profile_CIs_multiple.md)
  : Multi-value conditional GEV profile likelihood confidence intervals
- [`GEV_profile_loglik()`](https://opasche.github.io/ExtremeCI/reference/GEV_profile_loglik.md)
  : GEV profile log-likelihood
- [`GEV_maxlik()`](https://opasche.github.io/ExtremeCI/reference/GEV_maxlik.md)
  : Maximum-likelihood GEV estimate
- [`GEV_return_level()`](https://opasche.github.io/ExtremeCI/reference/GEV_return_level.md)
  : Compute return level from GEV parameters
- [`GEV_endpoint()`](https://opasche.github.io/ExtremeCI/reference/GEV_endpoint.md)
  : GEV endpoint
- [`GEV_log_likelihood()`](https://opasche.github.io/ExtremeCI/reference/GEV_log_likelihood.md)
  : GEV log likelihood
- [`GEV_change_parametrization()`](https://opasche.github.io/ExtremeCI/reference/GEV_change_parametrization.md)
  : GEV parameter vector reparametrization

## Plotting helpers

Helper functions for plotting the profile log-likelihood curves and
nonstationary quantile estimates with confidence intervals, as a
function of data, using `ggplot2`.

- [`plot_profile_loglik_curve()`](https://opasche.github.io/ExtremeCI/reference/plot_profile_loglik_curve.md)
  : Profile likelihood curve plot
- [`plot_data_quantile_ci()`](https://opasche.github.io/ExtremeCI/reference/plot_data_quantile_ci.md)
  : Non-stationary confidence bands plot
