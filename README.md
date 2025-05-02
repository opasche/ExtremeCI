
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExtremeCI

<!-- badges: start -->

<!-- badges: end -->

The `ExtremeCI` R package provides versatile algorithms to efficently
infer confidence intervals for extreme value statistics, using extreme
value theory extrapolation and the profile likelihood. These intervals
capture the uncertainty spread of extreme estimates more realistically
than most alternative methods, especially for return levels (high
quantiles) and for the shape parameter, whose uncertainties are
typically highly asymmetric (Coles, 2001,
<doi:10.1007/978-1-4471-3675-0>). For return levels and extreme
quantiles, the profile likelihood method requires reparametrization of
the likelihood function. With nonstationary models, this
reparametrization is nontrivial and requires repetition for each local
interval, which was not possible with alternative software. This package
provides a framework to construct both stationary and nonstationary
models with novel confidence endpoint search procedures based on binary
search, which do no require a prespecified range. The CIs can also be
inferred from weighted samples (work in progress). This package is
motivated by Zeder et al. (2023) <doi:10.1029/2023GL104090> and by
Pasche et al. (2026) <doi:10.48550/arXiv.2505.08578>.

## Installation

To install the development version of ExtremeCI from R, run

``` r
# install.packages("devtools")
devtools::install_github("opasche/ExtremeCI")
```

------------------------------------------------------------------------

Package created by Olivier C. PASCHE  
Research Institute for Statistics and Information Science,  
University of Geneva (CH), 2025.  
Supported by the Swiss National Science Foundation.
