% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{aspatial_siri}
\alias{aspatial_siri}
\title{Helper Function for Seasonal SIRI Simulation}
\usage{
aspatial_siri(
  initial_pop,
  season_length,
  mortality,
  transmission,
  recovery,
  fecundity,
  abundance_threshold,
  carrying_capacity,
  season
)
}
\arguments{
\item{initial_pop}{A vector of length 8 showing the initial abundance for
each combination of stage and compartment.}

\item{season_length}{The length of the season in days.}

\item{mortality}{A vector of length 8 with the mortality rates for each
stage and compartment in the season in question.}

\item{transmission}{A vector of length 8 with the transmission rates for each
stage in the season in question.}

\item{recovery}{A vector of length 8 with the recovery rates for each
infected stage in the season in question.}

\item{fecundity}{A vector of length 8 with the fecundity for each
reproductive segment.}

\item{abundance_threshold}{A quasi-extinction threshold below which a
population becomes extinct.}

\item{carrying_capacity}{A single numeric that indicates the carrying
capacity of the population in this season.}

\item{season}{Either "breeding" or "non-breeding."}
}
\value{
A vector of length 8 showing the abundance for each combination of
stage and compartment at the end of the season.
}
\description{
This function is an internal one that does the aspatial simulations within
one population for one timestep, for any given season.
}
\examples{
aspatial_siri(
 initial_pop = c(50000, 50000, 0, 1, 0, 0, 0, 0),
 season_length = 100,
 mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04),
 fecundity = c(0, 15/182, 0, 15/182, 0, 15/182, 0, 15/182),
 transmission = c(0.00002, 0.00001, 0, 0, 7.84e-06, 3.92e-06, 0, 0),
 recovery = c(0, 0, 0.05714286, 0.05714286, 0, 0, 0.1, 0.1),
 carrying_capacity = 150000,
 abundance_threshold = 10,
 season = "breeding"
)
}
