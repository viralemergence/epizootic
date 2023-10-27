#' Nested functions for a user-defined transformation of a population affected
#' by a disease outbreak.
#'
#' Modular functions for the disease simulator for performing a transformation
#' of a population across stages and disease compartments (and optionally
#' carrying capacity) at a specified time step via a user-defined function.
#'
#' @param replicates Number of replicate simulation runs.
#' @param time_steps Number of simulation time steps.
#' @param years_per_step Number of years per time step.
#' @param populations Number of populations.
#' @param demographic_stochasticity Boolean for optionally choosing demographic
#' stochasticity for the transformation.
#' @param density_stages Array of booleans or numeric (0,1) for each stage to
#' indicate which stages are affected by density.
#' @param abundance_threshold A quasi-extinction threshold below which a
#'  population becomes extinct.
#' @param mortality A vector of mortality rates, one for each
#'  combination of stages and compartments. Assumed by default to be daily
#'  mortality rates unless indicated otherwise (see below).
#' @param mortality_unit A vector indicating whether mortality rates are
#'  daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all 0.
#'  A list of vectors may be provided if this varies by season.
#' @param fecundity A vector of fecundity rates, one for each
#'  combination of stages and compartments for which fecundity applies (see
#'  \code{fecundity_mask} below).
#' @param fecundity_unit A vector indicating whether mortality rates are
#'  daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all 0.
#'  A list of vectors may be provided if this varies by season.
#' @param fecundity_mask A vector indicating which stages and
#'  compartments reproduce.
#' @param transmission A vector of transmission rates, one for each
#'  combination of stages and compartment for which transmission applies (see
#'  \code{transmission_mask} below.
#' @param transmission_unit A vector indicating whether mortality rates
#'  are daily or seasonal. 1 indicates seasonal, 0 indicates daily.
#' @param transmission_mask
#' @param transformation A user-defined function (optionally nested in a list with additional attributes) for performing transformation: \code{function(params)}, where \emph{params} is a list passed to the function containing:
#'   \describe{
#'     \item{\code{replicates}}{Number of replicate simulation runs.}
#'     \item{\code{time_steps}}{Number of simulation time steps.}
#'     \item{\code{years_per_step}}{Number of years per time step.}
#'     \item{\code{populations}}{Number of populations.}
#'     \item{\code{stages}}{Number of life cycle stages.}
#'     \item{\code{demographic_stochasticity}}{Boolean for optionally choosing demographic stochasticity for the transformation.}
#'     \item{\code{density_stages}}{Array of booleans or numeric (0,1) for each stage to indicate which stages are affected by density.}
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for each population at time step.}
#'     \item{\code{breeding_season_length}}{Array of breeding season lengths in days for each population at time step.}
#'     \item{\code{stage_abundance}}{Matrix of (current) abundance for each stage (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied at (current) time step.}
#'     \item{\code{simulator}}{\code{\link{SimulatorReference}} object with dynamically accessible \emph{attached} and \emph{results} lists.}
#'     \item{\code{additional attributes}}{Additional attributes when the transformation is optionally nested in a list.}
#'   }
#'   returns a transformed stage abundance matrix (or a list with stage abundance and carrying capacity)
#' @param simulator \code{\link{SimulatorReference}} object with dynamically accessible \emph{attached} and \emph{results} lists.
#' @param name Optional name for the transformation (default is "transformation").
#' @return Abundance (and capacity) transformation function: \code{function(r, tm, carrying_capacity, stage_abundance, occupied_indices)}, where:
#'   \describe{
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for each population at time step.}
#'     \item{\code{stage_abundance}}{Matrix of abundance for each stage (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied at time step.}
#'     \item{\code{returns}}{List with transformed stage abundance matrix (and optionally carrying capacity).}
#'   }
#' @export population_transformation
