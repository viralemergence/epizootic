
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
an *object-oriented* framework inside of R. R6 classes such as
`DiseaseModel` and `SimulationHandler` are used to store model
attributes, check them for consistency, pass them to parallel sessions
for simulation, and gather results and errors.

## Example

The core simulation engine of `epizootic` is the function
`disease_simulator`, which simulates spatially explicit disease dynamics
in populations. Here is the initial state of an idealized theoretical
scenario:

``` r
library(poems)
library(epizootic)
example_region <- Region$new(coordinates = data.frame(x = rep(seq(177.01, 177.05, 0.01), 5),
                             y = rep(seq(-18.01, -18.05, -0.01), each = 5)))
initial_abundance <- c(c(5000, 5000, 0, 1, 0, 0, 0, 0),
                          rep(c(5000, 5000, 0, 0, 0, 0, 0, 0), 24)) |>
      matrix(nrow = 8)
example_region$raster_from_values(initial_abundance[2,]) |>
  raster::plot(main = "Susceptible Adults")
```

<img src="man/figures/README-initial state-1.png" width="100%" />
