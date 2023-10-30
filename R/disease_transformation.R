#' Nested functions for a user-defined transformation of a population affected
#' by a disease outbreak.
#'
#' Modular functions for the disease simulator for performing a transformation
#' of a population across stages and disease compartments (and optionally
#' carrying capacity) at a specified time step via a user-defined function.
#'
#' @param replicates Number of replicate simulation runs.
#' @param time_steps Number of simulation time steps.
#' @param seasons Number of seasons.
#' @param populations Number of populations.
#' @param compartments \item{\code{compartments}}{Number of disease compartments
#'  (e.g., 3 for a SIR model). Default: 1.}
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
#' @param transmission_mask A vector indicating which stages and
#'  compartments are subject to transmission (i.e., classes susceptible to
#'  infection.) Must be the same length as \code{compartments}.
#' @param recovery A vector of recovery rates, one for each
#'  combination of stages and compartment for which recovery applies (see
#'  \code{recovery_mask} below.) If recovery varies by season, a list of
#'  recovery vectors the same length as \code{seasons} may be provided instead.
#' @param recovery_unit A vector
#'  indicating whether mortality rates are daily or seasonal. 1 indicates
#'  seasonal, 0 indicates daily.
#' @param recovery_mask A vector indicating which compartments are subject to 
#' recovery (i.e., infected classes that can recover.) Must be the same length as
#'  \code{compartments}.
#' @param transformation A user-defined function (optionally nested in a list
#' with additional attributes) for performing transformation:
#' \code{function(params)}, where \emph{params} is a list passed to the function
#' containing:
#'   \describe{
#'     \item{\code{replicates}}{Number of replicate simulation runs.}
#'     \item{\code{time_steps}}{Number of simulation time steps.}
#'     \item{\code{years_per_step}}{Number of years per time step.}
#'     \item{\code{populations}}{Number of populations.}
#'     \item{\code{stages}}{Number of life cycle stages.}
#'     \item{\code{compartments}}{Number of disease compartments (e.g., 3 for a
#'     SIR model).}
#'     \item{\code{demographic_stochasticity}}{Boolean for optionally choosing
#'     demographic stochasticity for the transformation.}
#'     \item{\code{density_stages}}{Array of booleans or numeric (0,1) for each
#'     stage to indicate which stages are affected by density.}
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
#'   returns a transformed stage abundance matrix (or a list with stage
#'   abundance and carrying capacity)
#' @param simulator \code{\link{SimulatorReference}} object with dynamically
#' accessible \emph{attached} and \emph{results} lists.
#' @param name Optional name for the transformation (default is
#' "transformation").
#' @return Abundance (and capacity) transformation function:
#' \code{function(r, tm, carrying_capacity, segment_abundance,
#' occupied_indices)}, where:
#'   \describe{
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'     each population at time step.}
#'     \item{\code{segment_abundance}}{Matrix of abundance for each
#'     stage-compartment combo (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied
#'     at time step.}
#'     \item{\code{returns}}{List with transformed stage abundance matrix (and
#'     optionally carrying capacity).}
#'   }
#' @export disease_transformation

disease_transformation <- function(replicates,
                                   time_steps,
                                   seasons,
                                   populations,
                                   demographic_stochasticity,
                                   density_stages,
                                   abundance_threshold,
                                   mortality,
                                   mortality_unit,
                                   fecundity,
                                   fecundity_unit,
                                   fecundity_mask,
                                   transmission,
                                   transmission_unit,
                                   transmission_mask,
                                   recovery,
                                   recovery_unit,
                                   recovery_mask,
                                   transformation,
                                   simulator,
                                   name = "transformation") {
  if (!is.null(transformation)) {
    # Derive stages
    stages <- length(density_stages)

    # Unpack transformation function and additional attributes from a list
    additional_attributes <- list()
    if (is.list(transformation)) {
      function_index <- which(unlist(lapply(transformation, is.function)))
      additional_attributes <- transformation[-function_index]
      transformation <- transformation[[function_index]]
    }

    params <- mget(ls())

    if (is.function(transformation)) {

      user_defined_function <- function(r, tm, carrying_capacity,
                                        segment_abundance, occupied_indices) {
        # Add attributes to be made available to the user-defined function
        params[["r"]] <- r
        params[["tm"]] <- tm
        params[["carrying_capacity"]] <- carrying_capacity
        params[["segment_abundance"]] <- segment_abundance
        params[["occupied_indices"]] <- occupied_indices

        # Run user-defined transformation function
        tryCatch({
          transformed <- transformation(params)
        },
        error = function(e){
          cli_abort(c("Error produced within user-defined {name} function:",
                      "x" = "{as.character(e)}"))
        })

        # Resolve and check returned transformed values
        if (is.list(transformed)) {
          if (!all(c("segment_abundance", "carrying_capacity") %in%
                   names(transformed))) {
            cli_abort(c("When returning a list, the user-defined {name} function
                      should contain segment_abundance and carrying_capacity.",
                      "x" = "List output contains {names(transformed)}."))
          }
        } else {
          # assume stage abundance matrix and place a list
          transformed <- list(segment_abundance = transformed)
        }

        # Warn if any negative or non-finite
        if (!all(is.finite(transformed[["segment_abundance"]]))) {
          cli_warn("Non-finite segment abundances returned by user-defined
                     {name} function at indices
                     {which(!is.finite(transformed[['segment_abundance']]))}.")
        }
        if ("carrying_capacity" %in% names(transformed) &&
            !all(is.finite(transformed[['carrying_capacity']]))) {
          cli_warn("Non-finite carrying capacities returned by user-defined
                     {name} function at indices
                     {which(!is.finite(transformed[['carrying_capacity']]))}.")
        }
        if (any(transformed[["segment_abundance"]][which(
          is.finite(transformed[['segment_abundance']])
          )] < 0)) {
          cli_warn(
            c(
              "Negative segment abundances returned by user-defined
                     {name} function at indices
                     {which(is.finite(transformed[['segment_abundance']] < 0))}"
            )
          )
        }
        if ("carrying_capacity" %in% names(transformed) &&
            any(transformed[["carrying_capacity"]][which(
              is.finite(transformed[['carrying_capacity']])
              )] < 0)) {
          cli_warn("Negative carrying capacities returned by user-defined
                     {name} function at indices
                     {which(is.finite(transformed[['carrying_capacity']] < 0))}"
          )
        }

        return(transformed)
      }
      environment(user_defined_function)[["name"]] <- name

      return(user_defined_function)
    }
  }
}

