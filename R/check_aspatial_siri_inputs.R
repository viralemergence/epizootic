#' Helper function to check the validity of inputs to the `siri_model_summer`
#' and `siri_model_winter` functions.
#'
#' This is an internal function that checks inputs to the `disease_simulator`
#' function to make sure they are valid, and sets default values for needed
#' inputs if their values are not supplied. The possible inputs for this
#' function are the same as the possible inputs to the `disease_simulator`
#' function.
#'
#' @param inputs A nested list with named elements:
#'   \describe{
#'     \item{\code{replicates}}{Number of replicate simulation runs.}
#'     \item{\code{time_steps}}{Number of simulation time steps.}
#'     \item{\code{populations}}{Number of populations.}
#'     \item{\code{stages}}{Number of life cycle stages.}
#'     \item{\code{compartments}}{Number of disease compartments (e.g., 3 for a
#'     SIR model).}
#'     \item{\code{abundance_threshold}}{A quasi-extinction threshold below
#'     which a population becomes extinct.}
#'     \item{\code{mortality}}{A vector of mortality rates, one for each
#'     combination of stages and compartments.}
#'     \item{\code{mortality_unit}}{A vector indicating whether mortality rates
#'     are daily or seasonal. 1 indicates seasonal, 0 indicates daily.}
#'     \item{\code{fecundity}}{A vector of fecundity rates, one for each
#'     combination of stages and compartments for which fecundity applies.}
#'     \item{\code{fecundity_unit}}{A vector indicating whether fecundity rates
#'     are daily or seasonal. 1 indicates seasonal, 0 indicates daily.}
#'     \item{\code{fecundity_mask}}{A vector indicating which stages and
#'     compartments reproduce.}
#'     \item{\code{transmission}}{A vector of transmission rates, one for each
#'     combination of stages and compartment for which transmission applies (see
#'     \code{transmission_mask} below.}
#'     \item{\code{transmission_unit}}{A vector indicating whether mortality
#'     rates are daily or seasonal. 1 indicates seasonal, 0 indicates daily.}
#'     \item{\code{transmission_mask}}{A vector indicating which stages and
#'     compartments are subject to transmission (i.e., classes susceptible to
#'     infection.)}
#'     \item{\code{recovery}}{A vector of recovery rates, one for each
#'     combination of stages and compartment for which recovery applies (see
#'     \code{recovery_mask} below.)}
#'     \item{\code{recovery_unit}}{A vector indicating whether mortality rates
#'     are daily or seasonal. 1 indicates seasonal, 0 indicates daily.}
#'     \item{\code{recovery_mask}}{A vector indicating which compartments are
#'     subject to recovery (i.e., infected classes that can recover.)}
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'     each population at time step.}
#'     \item{\code{breeding_season_length}}{Array of breeding season lengths in
#'     days for each population at time step.}
#'     \item{\code{segment_abundance}}{Matrix of (current) abundance for each
#'     stage-compartment combo (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied
#'     at (current) time step.}
#'     \item{\code{simulator}}{\code{\link{SimulatorReference}} object with
#'     dynamically accessible \emph{attached} and \emph{results} lists.}
#'     \item{\code{additional attributes}}{Additional attributes when the
#'     transformation is optionally nested in a list.}
#'   }
#'@return A list identical to the inputs (if there are no errors.)

check_aspatial_siri_inputs <- function(inputs) {
  list2env(inputs, environment())
  if (!is.null(inputs[["tm"]]) && !is.null(inputs[["time_steps"]])) {
    if (tm > time_steps) {
    cli_abort("You have indicated that we are at timestep {tm} even though
              total time steps is set at {time_steps}.")
    }
  }
  if (stages != 2) {
    cli_abort(c("The seasonal SIRI functions are built for a two-stage system.",
              "x" = "You have entered {stages} stage{?s}."))
  }
  if (compartments != 4) {
    cli_abort(c("The seasonal SIRI functions are built for 4 disease
                compartments.",
                "x" = "You have entered {compartments} compartment{?s}."))
  }
  if (is.null(fecundity)) {
    inputs[["fecundity"]] <- 0
  }
  if (is.null(fecundity_unit)) {
    inputs[["fecundity_unit"]] <- 0
  }
  if (is.null(fecundity_mask)) {
    inputs[["fecundity_mask"]] <- rep(0, 8)
  }
  if (!all(c(mortality[!as.logical(mortality_unit)],
             fecundity[!as.logical(fecundity_unit)],
             transmission[!as.logical(transmission_unit)],
             recovery[!as.logical(recovery_unit)]) <= 1) |
      !all(c(mortality, fecundity, transmission, recovery) >= 0)) {
    cli_abort("Daily mortality, fecundity, transmission, and recovery rates must
              be between 0 and 1 inclusive.")
  }
  if (sum(fecundity_mask) > 4) {
    cli_abort(c("The seasonal SIRI functions assume that only adults
                reproduce.",
              "x" = "Your fecundity mask indicates that juveniles reproduce."))
  }
  if (sum(transmission_mask) > 4) {
    cli_abort(c("The seasonal SIRI functions assume that only susceptibles and
                recovereds can be infected.",
                "x" = "Your transmission mask indicates that infecteds can be
                infected."))
  }
  if (sum(recovery_mask) > 4) {
    cli_abort(c("The seasonal SIRI functions assume that only two compartments
                (infected 1, infected 2+) can recover.",
                "x" = "Your recovery mask indicates that more than two
                compartments can recover."))
  }

  return(inputs)

}
