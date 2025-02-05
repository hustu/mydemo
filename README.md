
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mypkg <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->
<!-- badges: end -->

The goal of mydemo is to implement high dimensional mediation analysis
with confounding adjustment.

## Installation

You can install the development version of mydemo like so:

``` r
devtools::install_github("hustu/mydemo")
```

## Example

``` r
library(mydemo)
dat <- sim_data(n = 300, p = 2000)
X <- dat$X
Y <- dat$Y
M <- dat$M
COV <- dat$COV
```
