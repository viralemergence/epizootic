#' Study region for house finch conjunctivitis simulations
#'
#' This is a raster that defines the spatial extent and resolution of the house
#' finch conjunctivitis example in the vignette. It encompasses terrestrial
#' North America from southern Mexico to southern Canada. It has a resolution of
#'  46.375 by 46.375 kilometers, adding up to 6,355 possible populations of
#'  house finches. The projected coordinate system for the study region is
#'  Albers equal-area conic, with a reference longitude of 94.5 W and standard
#'  parallels at 21.5 N and 47.5 N.
#' @format A raster object with 1 layer and 17066 cells defining the spatial
#' scope of a simulation.
"finch_region"

#' Raster of breeding season length for the house finch
#'
#' This is a RasterBrick object containing data on the breeding season
#' length in days of the house finch in North America from 1994 to 2016. I
#' created this dataset using a machine learning algorithm on season length data
#' from the eastern bluebird, which is noted to have a very similar breeding
#' season to the house finch. The raster has the same resolution as the
#' \code{\link{finch_region}}.
#' @format A RasterBrick with 17066 cells and 23 layers.
"bsl_raster"

#' Initial house finch abundance
#'
#' This numeric matrix contains initial abundances of house finches in 1994 for
#' the house finch conjunctivitis vignette. These abundances were themselves
#' simulated using `epizootic` and do not represent empirical estimates of house
#' finch abundance in 1994. The matrix has 6,355 columns, one for each
#' population in the [finch_region], and 8 rows, one for each combination of
#' life cycle stage and disease compartment (row 1: susceptible juveniles, row
#' 2: susceptible adults, row 3: juveniles infected for the first time, row 4:
#' adults infected for the first time, row 5: recovered juveniles, row 6:
#' recovered adults, row 7: re-infected juveniles, row 8: re-infected adults.)
#' This initial abundance matrix has only zeroes in rows 3-8 because the disease
#' has not yet broken out at the start of the simulation.
#' @format A numeric matrix with 8 rows and 6355 columns.
"initial_abundance"

#' House finch habitat suitability
#'
#' This is a RasterStack containing data on habitat suitability for
#' the house finch in North America from 1994 to 2016. This habitat suitability
#' stack was generated using a species distribution model. The predictors for
#' the SDM were the 12 bioclimatic variables, plus an urbanization index
#' (proportion of each grid cell with urban land use.) The occurrences for the
#' SDM came from quality-checked GBIF records.
#'
#' The raster has the same resolution and projection as the [finch_region].
#' Habitat suitability ranges from 0 to 1, with 1 being the most suitable.
#' @format A RasterStack with 17066 cells and 23 layers.
"habitat_suitability"
