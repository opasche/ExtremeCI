# ExtremeCI

The `ExtremeCI` R package provides versatile algorithms to efficiently
infer confidence intervals for extreme value statistics, such as extreme
quantiles and return levels, that are representative of the asymmetric
uncertainty spread, using extreme value theory extrapolation and the
profile likelihood (see e.g., [Coles,
2001](https://doi.org/10.1007/978-1-4471-3675-0)). Unlike existing
algorithms, the CI endpoints are found without the need for a strict
prespecified range, can be covariate-dependent, and can be based on
weighted samples. This package is motivated by [Zeder et
al. (2023)](https://doi.org/10.1029/2023GL104090) and by [Pasche et
al. (2026)](https://doi.org/10.1007/s10687-026-00536-9).

## Installation

To install the latest release of `ExtremeCI` from CRAN, in an R session,
simply run:

``` r

install.packages("ExtremeCI")
```

Or, to install the current development version of `ExtremeCI`, in an R
session, run:

``` r

# install.packages("devtools")
devtools::install_github("opasche/ExtremeCI")
```

## Motivation

Confidence intervals based on the profile-likelihood method capture the
uncertainty spread of extreme value statistics more realistically than
most alternative methods (e.g. based on bootstrapping or the
delta-method), especially for return levels (high quantiles) and for the
shape parameter, whose uncertainties are typically highly asymmetric
(see e.g., [Coles, 2001](https://doi.org/10.1007/978-1-4471-3675-0)).
For return levels and extreme quantiles, the profile likelihood method
requires reparametrization of the likelihood function. With
nonstationary models, this reparametrization is nontrivial and requires
repetition for each local interval, which was not possible with
alternative software. Furthermore, preexisting implementations require
prespecified ranges for the CI endpoints search, which is impractical
when the uncertainty spread is large or covariate-dependent, or when
used as step in a greater pipeline. The `ExtremeCI` R package provides a
framework to construct both stationary and nonstationary models with
binary-search procedures to efficiently find the profile likelihood CI
endpoints, without requiring a prespecified range. Lastly, weighted
inference can optionally be performed using sample weights. In summary,
the novelties of this package for realistic extreme value confidence
intervals are: - efficient profile-likelihood CI endpoint search
procedures without prespecified ranges, for both the peaks-over
threshold (GPD) and block maxima (GEV) approaches; - covariate-dependent
CIs for nonstationary extreme value models, including for extreme
quantiles and return levels; - the option to use sample weights in
inference.

## References

Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
Values*. Springer.
<doi:%5B10.1007/978-1-4471-3675-0>\](<https://doi.org/10.1007/978-1-4471-3675-0>).

Pasche, O. C., Lam, H., and Engelke, S. (2026). “Extreme Conformal
Prediction: Reliable Intervals for High-Impact Events.” *Extremes*.
<doi:%5B10.1007/s10687-026-00536-9>\](<https://doi.org/10.1007/s10687-026-00536-9>).

Zeder, J., Sippel, S., Pasche, O. C., Engelke, S., and Fischer, E. M.
(2023). “The effect of a short observational record on the statistics of
temperature extremes.” *Geophysical Research Letters* 50(16),
e2023GL104090.
<doi:%5B10.1029/2023GL104090>\](<https://doi.org/10.1029/2023GL104090>).

------------------------------------------------------------------------

Package created by Olivier C. PASCHE  
Research Institute for Statistics and Information Science,  
University of Geneva (CH), 2025.  
Supported by the Swiss National Science Foundation.
