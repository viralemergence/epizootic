
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epizootic

<!-- badges: start -->

[![R-CMD-check](https://github.com/viralemergence/epizootic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/viralemergence/epizootic/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`epizootic` is an extension to `poems`, a spatially-explicit,
process-explicit, pattern-oriented framework for modeling population
dynamics. This extension adds functionality for modeling disease
dynamics in wildlife. It also adds capability for seasonality and for
unique dispersal dynamics for each life cycle stage.

## Installation

You can install the development version of epizootic from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("viralemergence/epizootic")
```

## About R6 classes

`poems` and `epizootic` run on
[R6](https://r6.r-lib.org/articles/Introduction.html) classes. R is
primarily a *functional* programming language in which the primary units
of programming are expressions and functions. Here we use R6 to create
an *object-oriented* framework inside of R.
