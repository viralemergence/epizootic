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
#'     \item{\code{abundance_threshold}}{A quasi-extinction threshold at
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
#'     \item{\code{simulator}}{[`poems::SimulatorReference`] object with
#'     dynamically accessible \emph{attached} and \emph{results} lists.}
#'     \item{\code{additional attributes}}{Additional attributes when the
#'     transformation is optionally nested in a list.}
#'   }
#'@return A list identical to the inputs (if there are no errors.)

check_aspatial_siri_inputs <- function(inputs) {
  list2env(inputs, environment())

  segments <- inputs[["stages"]] * inputs[["compartments"]]
  if (any(segment_abundance < 0)) {
    cli_abort("x" = "Abundance values must be non-negative.",
              "i" = "Populations {which(apply(segment_abundance, 2, function(x) any(x < 0)))}
                     contain negative values at time step {tm}.")
  }

  if (any(is.na(segment_abundance))) {
    cli_abort("x" = "Abundance values must be non-missing.",
              "i" = "Populations {which(apply(is.na(segment_abundance), 2, function(x) any(x)))}
                     contain missing values at time step {tm}.")
  }

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
  if (is.null(inputs[["fecundity"]])) {
    inputs[["fecundity"]] <- fecundity <- 0
  }
  if (is.null(inputs[["fecundity_unit"]])) {
    inputs[["fecundity_unit"]] <- fecundity_unit <- 0
  }
  if (is.null(inputs[["fecundity_mask"]])) {
    inputs[["fecundity_mask"]] <- fecundity_mask <- rep(0, 8)
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

  # Check fecundity and survival
  mortality <- inputs[["mortality"]]
  if (is.vector(mortality)) {
    if (length(mortality) != segments) {
      cli_abort(
        c("The {.var mortality} vector must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment.")
      )
    }
  } else {
    cli_abort(
      c("{.var mortality} must be a vector, not
        {.obj_type_friendly {mortality}}.")
    )
  }
  if (is.null(inputs[["mortality_unit"]])) {
    inputs[["mortality_unit"]] <- rep(0, segments)
  } else {
    unit_values <- inputs[["mortality_unit"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("mortality_unit values must be 0 or 1"))
    }
    if (length(unit_values) == 1) {
      inputs[["mortality_unit"]] <- rep(unit_values, 8)
    }
    if (length(mortality) != length(inputs[["mortality_unit"]])) {
      mortality_unit <- inputs[["mortality_unit"]]
      cli_abort(c("`mortality` and `mortality_unit` must be the same length.",
                  "*" = "`mortality` is length {length(mortality)}.",
                  "*" = "`mortality_unit` is length {length(mortality_unit)}."))
    }
  }

  apply_mask <- function(transmission_unit, transmission_mask) {
    if (is.vector(transmission_unit) && is.vector(transmission_mask)) {
      unit_index <- 1
      result <- integer(length(transmission_mask))
      for (i in seq_along(transmission_mask)) {
        if (transmission_mask[i] == 1) {
          result[i] <- transmission_unit[unit_index]
          unit_index <- unit_index + 1
        } else {
          result[i] <- 0
        }
      }
      result
    } else {
      cli_abort(c("Please only input vectors for transmission_unit and
                   transmission_mask."))
    }
  }

  fecundity <- inputs[["fecundity"]]
  if (is.vector(fecundity)) {
    if (length(fecundity) != segments) {
      if (is.vector(inputs[["fecundity_mask"]])) {
        fecundity_mask <- inputs[["fecundity_mask"]]
        z <- rep(0, segments)
        multiplier <- sum(fecundity_mask) / length(fecundity)
        z[as.logical(fecundity_mask)] <- rep(fecundity, multiplier)
        fecundity <- inputs[["fecundity"]] <- z
      } else {
        cli_abort(
          c("The {.var fecundity} vector must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment. Or, provide a fecundity mask so
          that fecundity values may be assigned to the appropriate stages and
          compartments.")
        )
      }
    }
  } else {
    cli_abort(
      c("{.var fecundity} must be a vector, not
        {.obj_type_friendly {fecundity}}.")
    )
  }

  if (is.null(inputs[["fecundity_mask"]])) {
    inputs[["fecundity_mask"]] <- rep(1, segments)
  } else {
    unit_values <- inputs[["fecundity_mask"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("fecundity_mask values must be 0 or 1"))
    }
  }
  fecundity <- inputs[["fecundity"]] <- map2_dbl(fecundity,
                                             inputs[["fecundity_mask"]],
                                             ~ ifelse(.y == 0, 0, .x))

  if (is.null(inputs[["fecundity_unit"]])) {
    inputs[["fecundity_unit"]] <- rep(0, segments)
  } else {
    unit_values <- inputs[["fecundity_unit"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("fecundity_unit values must be 0 or 1. Unit values are
                  {unit_values}."))
    }
    if (length(fecundity) != length(inputs[["fecundity_unit"]])) {
      
      if (length(inputs[["fecundity_unit"]]) == 1) {
        inputs[["fecundity_unit"]] <- rep(inputs[["fecundity_unit"]], segments)
      } else if (length(inputs[["fecundity_unit"]]) == sum(inputs[["fecundity_mask"]])) {
        inputs[["fecundity_unit"]] <- apply_mask(inputs[["fecundity_unit"]], 
                                                inputs[["fecundity_mask"]])
      } else {
        fecundity_unit <- inputs[["fecundity_unit"]]
        cli_abort(c("`fecundity` and `fecundity_unit` must be the same length.",
                    "*" = "`fecundity` is length {length(fecundity)}.",
                    "*" = "`fecundity_unit` is length {length(fecundity_unit)}."))
      }
    }
  }

  # Check transmission and recovery
  transmission <- inputs[["transmission"]]
  if (is.vector(transmission)) {
    if (length(transmission) != segments && !is.vector(inputs[["transmission_mask"]])) {
      cli_abort(
        c("The {.var transmission} vector must have
        {inputs[['stages']]*inputs[['compartments']]} elements, one for each
        combination of stage and compartment. Or, provide a transmission mask so
        that transmission values may be assigned to the appropriate stages and
        compartments.")
      )
    }
  } else {
    cli_abort(
      c("{.var transmission} must be a vector or list, not
        {.obj_type_friendly {transmission}}.")
    )
  }
  if (is.null(inputs[["transmission_unit"]])) {
    inputs[["transmission_unit"]] <- rep(0, segments)
  } else {
    unit_values <- inputs[["transmission_unit"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("transmission_unit values must be 0 or 1. Unit values are
                  {unit_values}."))
    }
    if (length(transmission) != length(inputs[["transmission_unit"]])) {
      transmission_unit <- inputs[["transmission_unit"]]
      cli_abort(c("`transmission` and `transmission_unit` must be the same length.",
                  "*" = "`transmission` is length {length(transmission)}.",
                  "*" = "`transmission_unit` is length {length(transmission_unit)}."))
    }
  }
  if (is.null(inputs[["transmission_mask"]])) {
    inputs[["transmission_mask"]] <- rep(0, segments) |> replace(1:stages, 1)
  } else {
    unit_values <- inputs[["transmission_mask"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("transmission_mask values must be 0 or 1"))
    }
  }
  if (length(transmission) != segments &&
      length(inputs[['transmission_mask']]) == segments) {
    z <- rep(0, length(inputs[['transmission_mask']]))
    z[as.logical(inputs[['transmission_mask']])] <- transmission
    transmission <- inputs[["transmission"]] <- z
  }
  # use the function to apply the mask to the transmission unit
  inputs[["transmission_unit"]] <- apply_mask(inputs[["transmission_unit"]],
                                              inputs[["transmission_mask"]])
  inputs[["transmission"]] <- map2_dbl(inputs[['transmission']],
                                   inputs[["transmission_mask"]],
                                   ~ ifelse(.y == 0, 0, .x))

  recovery <- inputs[["recovery"]]
  if (is.null(recovery)) {
    recovery <- inputs[["recovery"]] <- rep(0, segments)
  }
  if (is.vector(recovery)) {
    if (length(recovery) != segments) {
      if (is.vector(inputs[["recovery_mask"]])) {
        recovery_mask <- inputs[["recovery_mask"]]
        z <- rep(0, segments)
        multiplier <- sum(recovery_mask) / length(recovery)
        z[as.logical(recovery_mask)] <- rep(recovery, multiplier)
        recovery <- inputs[["recovery"]] <- z
      } else {
        cli_abort(
          c("The {.var recovery} vector must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment. Or, provide a recovery mask so
          that recovery values may be assigned to the appropriate stages and
          compartments.")
        )
      }
    }
  } else {
    cli_abort(
      c("{.var recovery} must be a vector, not
        {.obj_type_friendly {recovery}}.")
    )
  }
  if (is.null(inputs[["recovery_unit"]])) {
    inputs[["recovery_unit"]] <- rep(0, segments)
  } else {
    unit_values <- inputs[["recovery_unit"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("recovery_unit values must be 0 or 1. Unit values are
                  {unit_values}."))
    }
    if (length(recovery) != length(inputs[["recovery_unit"]])) {
      recovery_unit <- inputs[["recovery_unit"]]
      if (length(inputs[["recovery_unit"]]) == 1) {
        inputs[["recovery_unit"]] <- rep(inputs[["recovery_unit"]], segments)
      } else if (length(inputs[["recovery_unit"]]) == sum(inputs[["recovery_mask"]])) {
        inputs[["recovery_unit"]] <- apply_mask(inputs[["recovery_unit"]], 
                                                inputs[["recovery_mask"]])
      } else {
        cli_abort(c("`recovery` and `recovery_unit` must be the same length.",
                  "*" = "`recovery` is length {length(recovery)}.",
                  "*" = "`recovery_unit` is length {length(recovery_unit)}."))
      }
    }
  }
  if (is.null(inputs[["recovery_mask"]])) {
    if (inputs[["compartments"]] > 1) {
      recovery_mask <- rep(0, segments) |> replace((stages + 1):(stages * 2), 1)
    } else if (inputs[["compartments"]] == 1) {
      recovery_mask <- rep(0, segments)
    }
    inputs[["recovery_mask"]] <- recovery_mask
  } else {
    unit_values <- inputs[["recovery_mask"]] |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("recovery_mask values must be 0 or 1"))
    }
  }
  # use the function to apply the mask to the recovery unit
  inputs[["recovery_unit"]] <- apply_mask(inputs[["recovery_unit"]],
                                          inputs[["recovery_mask"]])
  recovery <- inputs[["recovery"]] <- map2_dbl(recovery,
                                           inputs[["recovery_mask"]],
                                           ~ ifelse(.y == 0, 0, .x))

  return(inputs)
}
