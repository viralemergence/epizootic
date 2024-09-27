#' Helper function to check the validity of inputs to the `disease_simulator`
#' function.
#'
#' This is an internal function that checks inputs to the `disease_simulator`
#' function to make sure they are valid, and sets default values for needed
#' inputs if their values are not supplied. The possible inputs for this
#' function are the same as the possible inputs to the `disease_simulator`
#' function.
#'
#'@param inputs Nested list/object with named elements:
#'@param inputs Nested list/object with named elements:
#'\describe{
#'  \item{\code{random_seed}}{Number to seed the random number generation for
#'  stochasticity.}
#'  \item{\code{replicates}}{Number of replicate simulation runs (default is
#'  1.)}
#'  \item{\code{time_steps}}{Number of simulation years. Required input.}
#'  \item{\code{seasons}}{Number of seasons per year (default is 2.)}
#'  \item{\code{populations}}{Number of populations. Required input.}
#'  \item{\code{coordinates}}{Data frame (or matrix) of X-Y population
#'  coordinates.}
#'  \item{\code{stages}}{Number of life cycle stages. Default: 1.}
#'  \item{\code{compartments}}{Number of disease compartments (e.g., 3 for a
#'  SIR model). Default: 1.}
#'  \item{\code{region}}{A [`poems::Region`] object
#'  defining the study region.}
#'  \item{\code{initial_abundance}}{Array (or matrix) of initial abundances.
#'  There must be one column per population and one row per compartment/stage
#'  combination. By default, this should be in the order compartment by stage,
#'  e.g., 2 stage classes plus a SI model should be ordered as S1, S2, I1, I2.
#'  If a region object is attached, then initial abundance may be provided in
#'  the form of a raster with the same specs as the region raster and one
#'  layer per stage/compartment combination. If there is only one
#'  stage/compartment combination you may provide a vector with length
#'  \code{populations}. Required input.}
#'  \item{\code{carrying_capacity}}{Array (matrix) of carrying capacity values
#'  at each population cell (\code{populations} rows by \code{time_steps}
#'  columns when across time). Required input.}
#'  \item{\code{breeding_season_length}}{Array (matrix) of breeding season
#'  length values in days at each population cell (\code{populations} rows by
#'  \code{time_steps} columns when across time). Can also be a vector of length
#'  \code{populations} if the breeding season length does not change over
#'  time.}
#'  \item{\code{season_lengths}}{Vector of season lengths in days. Length must
#'  equal \code{seasons}. If neither \code{breeding_season_length} nor
#'  \code{season_lengths} are provided, season lengths will default to
#'  \code{365/seasons}.}
#'  \item{\code{correlation}}{List containing either an environmental
#'  correlation matrix (correlation_matrix), a pre-calculated transposed
#'  (Cholesky) decomposition matrix (t_decomposition_matrix), or a compact
#'  transposed (Cholesky) decomposition matrix (t_decomposition_compact_matrix)
#'  and a corresponding map of population indices (t_decomposition_compact_map),
#'  as per [`poems::SpatialCorrelation`] class attributes.}
#'  \item{\code{mortality}}{A vector of mortality rates, one for each
#'  combination of stages and compartments. Assumed by default to be daily
#'  mortality rates unless indicated otherwise (see below). If mortality varies
#'  by season, a list of mortality vectors with the same length as \code{seasons}
#'  may be provided instead. Required input.}
#'  \item{\code{mortality_unit}}{A vector indicating whether mortality rates are
#'  daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all 0.
#'  A list of vectors may be provided if this varies by season.}
#'  \item{\code{fecundity}}{A vector of fecundity rates, one for each
#'  combination of stages and compartments for which fecundity applies (see
#'  \code{fecundity_mask} below). If fecundity varies among seasons, a list of
#'  fecundity vectors with the same length as \code{seasons} may be provided.
#'   Required input.}
#'  \item{\code{fecundity_unit}}{A vector indicating whether fecundity rates are
#'  daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all 0.
#'  A list of vectors may be provided if this varies by season.}
#'  \item{\code{fecundity_mask}}{A vector indicating which stages and
#'  compartments reproduce. Must be the same length as
#'  \code{stages * compartments}. A list of vectors may be provided if this
#'  varies by season. If no fecundity mask is provided, then it is assumed that
#'  all stages and compartments reproduce.}
#'  \item{\code{abundance_threshold}}{A quasi-extinction threshold at which a
#'  population becomes extinct. Default: 0.}
#'  \item{\code{demographic_stochasticity}}{Boolean for choosing demographic
#'  stochasticity for transition, dispersal, and/or other processes (default is
#'  TRUE).}
#'  \item{\code{transmission}}{A vector of transmission rates, one for each
#'  combination of stages and compartment for which transmission applies (see
#'  \code{transmission_mask} below). If transmission varies by season, a list of
#'  transmission vectors with the same length as \code{seasons} may be provided
#'  instead. Required input.}
#'  \item{\code{transmission_unit}}{A vector indicating whether transmission
#'  is daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all
#'  0. A list of vectors may be provided if this varies by season.}
#'  \item{\code{transmission_mask}}{A vector indicating which stages and
#'  compartments are subject to transmission (i.e., classes susceptible to
#'  infection.) Must be the same length as \code{compartments}. A list of
#'  vectors may be provided if this varies by season. If no transmission mask is
#'   provided, then it is assumed that all stages in the first compartment are
#'  susceptible to infection.}
#'  \item{\code{recovery}}{A vector of recovery rates, one for each
#'  combination of stages and compartment for which recovery applies (see
#'  \code{recovery_mask} below.) If recovery varies by season, a list of
#'  recovery vectors the same length as \code{seasons} may be provided instead.}
#'  \item{\code{recovery_unit}}{A vector
#'  indicating whether recovery rates are daily or seasonal. 1 indicates
#'  seasonal, 0 indicates daily. Default: all 0. A list of vectors may be
#'  provided if this varies by season.}
#'  \item{\code{recovery_mask}}{A
#'  vector indicating which compartments are subject to recovery (i.e., infected
#'  classes that can recover.) Must be the same length as \code{compartments}.
#'  A list of vectors may be provided if this varies by season. If no recovery
#'  mask is provided, then it is assumed that all stages in the second
#'  compartment can recover, if there is a second compartment.}
#'  \item{\code{dispersal}}{A list that is either length 1 or the same length as
#'  \code{stages}. If it is length 1, the same dispersal will be applied across
#'  all stages. Within each element of the list, there should be either a function,
#'  a matrix of dispersal rates between populations (source columns to target
#'  rows) or a list of data frames of non-zero dispersal rates and indices for
#'  constructing a compact dispersal matrix, and optional changing rates over
#'  time (as per class [`poems::DispersalGenerator`] \emph{dispersal_data}
#'  attribute).}
#'  \item{\code{dispersal_source_n_k}}{Dispersal proportion (p) density
#'  dependence via source population abundance divided by carrying capacity
#'  (n/k), where p is reduced via a linear slope (defined by two list items)
#'  from n/k <= \emph{cutoff} (p = 0) to n/k \>= \emph{threshold} (aliases:
#'  \emph{dispersal_n_k_cutoff} & \emph{dispersal_n_k_threshold}).}
#'  \item{\code{dispersal_target_k}}{Dispersal rate (r) density dependence via
#'  target population carrying capacity (k), where r is reduced via a linear
#'  slope (through the origin) when k <= \emph{threshold} (alias:
#'  \emph{dispersal_k_threshold}).}
#'  \item{\code{dispersal_target_n}}{Dispersal
#'  rate (r) density dependence via target population abundance (n), where r is
#'  reduced via a linear slope (defined by two list items) from n \>=
#'  \emph{threshold} to n <= \emph{cutoff} (r = 0) or vice versa (aliases:
#'  \emph{dispersal_n_threshold} & \emph{dispersal_n_cutoff}).}
#'  \item{\code{dispersal_target_n_k}}{Dispersal rate (r) density dependence via
#'  target population abundance divided by carrying capacity (n/k), where r is
#'  reduced via a linear slope (defined by two list items) from n/k \>=
#'  \emph{threshold} to n/k <= \emph{cutoff} (r = 0) or vice versa.}
#'  \item{\code{season_functions}}{A list of population transformation functions
#'  (functions that change abundance across stages and compartments) the same
#'  length as \code{seasons.} The function must be in the form
#'  \code{function(params)}, where \code{params} is a list passed to the
#'  function containing:
#'       \describe{
#'          \item{\code{replicates}}{Number of replicate simulation runs
#'          (default is 1.)}
#'          \item{\code{time_steps}}{Number of simulation years. Required
#'          input.}
#'          \item{\code{seasons}}{Number of seasons per year (default is 2.)}
#'          \item{\code{populations}}{Number of populations. Required input.}
#'          \item{\code{stages}}{Number of life cycle stages. Default: 1.}
#'          \item{\code{compartments}}{Number of disease compartments (e.g., 3
#'          for a SIR model). Default: 1.}
#'          \item{\code{breeding_season_length}}{Array (matrix) of breeding
#'          season length values in days at each population cell
#'          (\code{populations} rows by \code{time_steps} columns when across
#'          time).}
#'          \item{\code{season_lengths}}{Vector of season lengths in days.
#'          Length must equal \code{seasons}.}
#'          \item{\code{mortality}}{A vector of mortality rates, one for each
#'          combination of stages and compartments. Assumed by default to be
#'          daily mortality rates. Required input.}
#'          \item{\code{mortality_unit}}{A vector indicating whether mortality
#'          rates are daily or seasonal. 1 indicates seasonal, 0 indicates
#'          daily. Default: all 0.}
#'          \item{\code{fecundity}}{A vector of fecundity rates, one for each
#'          combination of stages and compartments for which fecundity applies
#'          (see \code{fecundity_mask} below). Required input.}
#'          \item{\code{fecundity_unit}}{A vector indicating whether mortality
#'          rates are daily or seasonal. 1 indicates seasonal, 0 indicates
#'          daily. Default: all 0.}
#'          \item{\code{fecundity_mask}}{A vector indicating which stages and
#'          compartments reproduce. Must be the same length as
#'          \code{stages * compartments}.}
#'          \item{\code{abundance_threshold}}{A quasi-extinction threshold below
#'          which a population becomes extinct. Default: 0.}
#'          \item{\code{demographic_stochasticity}}{Boolean for choosing
#'          demographic stochasticity for transition, dispersal, and/or other
#'          processes (default is TRUE).}
#'          \item{\code{transmission}}{A vector of transmission rates, one for
#'          each combination of stages and compartments. Assumed by default to
#'          be daily transmission rates. Required input.}
#'          \item{\code{transmission_unit}}{A vector indicating whether
#'          mortality rates are daily or seasonal. 1 indicates seasonal, 0
#'          indicates daily. Default: all 0.}
#'          \item{\code{recovery}}{A vector of recovery rates, one for each
#'          combination of stages and compartment for which recovery applies
#'          (see \code{recovery_mask} below.)}
#'          \item{\code{recovery_unit}}{A vector indicating whether mortality
#'          rates are daily or seasonal. 1 indicates seasonal, 0 indicates
#'          daily. Default: all 0.}
#'          \item{\code{recovery_mask}}{A vector indicating which compartments
#'          are subject to recovery (i.e., infected classes that can recover.)
#'          Must be the same length as \code{compartments}.}
#'          \item{\code{r}}{Simulation replicate.}
#'          \item{\code{tm}}{Simulation time step.}
#'          \item{\code{carrying_capacity}}{Array of carrying capacity values
#'          for each population at time step.}
#'          \item{\code{segment_abundance}}{Matrix of abundance for each
#'          combination of stage and compartment (rows) and population (columns)
#'           at time step.}
#'          \item{\code{occupied_indices}}{Array of indices for populations
#'          occupied at time step.}
#'          \item{\code{simulator}}{[`poems::SimulatorReference`] object
#'          with dynamically accessible \emph{attached} and \emph{results}
#'          lists.}
#'          \item{\code{additional attributes}}{Additional attributes when the
#'          transformation is optionally nested in a list.}
#'          }
#'       and returns a transformed stage abundance matrix.}
#'  \item{\code{simulation_order}}{A list the same length as \code{seasons}.
#'  Each element in the list is a vector of named simulation processes in the
#'  desired order. Processes must be one of "transition", "dispersal",
#'  "season_functions", or "results."}
#'  \item{\code{dispersal_type}}{A character vector that may contain "pooled"
#'  (if all individuals disperse the same), "stages", "compartments", or
#'  "segments", if different stages, compartments, or stage-compartment
#'  combinations disperse differently. If "pooled" is chosen,
#'  \code{dispersal} must be a list of length 1. If "stages" is chosen, it
#'  must be the same length as \code{stages}, if "compartments" is chosen,
#'  it must be the same length as \code{compartments}, and if "segments" is
#'  chosen, it must be the same length as stages*compartments. The default
#' value is "pooled".}
#'  \item{\code{results_selection}}{List of results selection from: "abundance"
#'  (default), "ema", "extirpation", "extinction_location",
#'  "occupancy"; "summarize" (default) or "replicate".}
#'  \item{\code{results_breakdown}}{A string with one of these values:
#' "segments" (default),
#' "compartments", "stages" or "pooled." "segments" returns results for each
#' segment (stage x compartment combination.) "compartments" returns results for
#' each disease compartment. "stages" returns results for each life cycle stage.
#' "pooled" returns results that are not broken down by stage or compartment.}
#' \item{\code{verbose}}{TRUE or FALSE, indicating if the user wants informative
#' messages throughout the simulation process.}
#'}
#'@return A list identical to the inputs, except with default values supplied
#' to fill in any crucial missing values, as explained in the documentation
#' above.
#'@export

check_simulator_inputs <- function(inputs) {

  # Unpack inputs if they are in the form of a DiseaseModel
  if (inherits(inputs, "DiseaseModel")) {
    if (!inputs$is_consistent()) {
      cli_abort(c("The `DiseaseModel` input contains the following inconsistent
                  attributes:",
                  "x" = "{inputs$inconsistent_attributes()}"))
    } else if (!inputs$is_complete()) {
      cli_abort(c("The `DiseaseModel` input contains the following incomplete
                  attributes:",
                  "x" = "{inputs$incomplete_attributes()}"))
    } else {
      inputs <- inputs$get_attributes()
    }
  }

  # Stop if minimal inputs are not present
  if (is.null(inputs[["time_steps"]]) || is.null(inputs[["populations"]]) ||
      is.null(inputs[["initial_abundance"]]) ||
      is.null(inputs[["carrying_capacity"]]) ||
      is.null(inputs[["mortality"]]) || is.null(inputs[["transmission"]]) ||
      is.null(inputs[["fecundity"]]) || is.null(inputs[["simulation_order"]])) {
    incomplete_inputs <- if (is.null(inputs[["time_steps"]]))
      "time_steps"
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["populations"]]))
        "populations"
    )
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["initial_abundance"]]))
        "initial_abundance"
    )
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["carrying_capacity"]]))
        "carrying_capacity"
    )
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["fecundity"]]))
        "fecundity"
    )
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["mortality"]]))
        "mortality"
    )
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["transmission"]]))
        "transmission"
    )
    incomplete_inputs <- c(
      incomplete_inputs, if (is.null(inputs[["simulation_order"]]))
        "simulation_order"
    )
    cli_abort(c("Minimal inputs required to run simulation are missing.",
                "x" = "Your input does not include {incomplete_inputs}."))
  }

  ## Unpack inputs and calculate/initialize re-usable variables ##

  # General simulation settings
  if (!is.null(inputs[["random_seed"]])) {
    set.seed(inputs[["random_seed"]])
  }
  inputs[["replicates"]] <- ifelse(
    is.null(inputs[["replicates"]]), 1, inputs[["replicates"]]
  )
  inputs[["seasons"]] <- ifelse(
    is.null(inputs[["seasons"]]), 2, inputs[["seasons"]]
  )
  inputs[["stages"]] <- ifelse(
    is.null(inputs[["stages"]]), 1, inputs[["stages"]]
  )
  stages <- inputs[["stages"]]
  inputs[["compartments"]] <- ifelse(
    is.null(inputs[["compartments"]]), 1, inputs[["compartments"]]
  )
  compartments <- inputs[["compartments"]]
  inputs[["demographic_stochasticity"]] <- ifelse(
    is.null(inputs[["demographic_stochasticity"]]),
    TRUE,
    inputs[["demographic_stochasticity"]]
  )

  # Initial abundance for each stage and population
  if (!is.null(inputs[["region"]]) &&
    inputs[["region"]][["use_raster"]] &&
    any(c("RasterLayer", "RasterStack", "RasterBrick") %in%
    class(inputs[["initial_abundance"]]))) {

      # Convert the RasterBrick to a matrix
      inputs[["initial_abundance"]] <- raster::values(inputs[["initial_abundance"]])

      # Subset the matrix
      inputs[["initial_abundance"]] <- inputs[["initial_abundance"]][inputs[["region"]][["region_indices"]],]
    } else if (is.vector(inputs[["initial_abundance"]]) &&
             is.numeric(inputs[["initial_abundance"]])) {
    inputs[["initial_abundance"]] <- matrix(inputs[["initial_abundance"]],
                                            nrow = 1)
  } else if (is.array(inputs[["initial_abundance"]])) {
    inputs[["initial_abundance"]] <- inputs[["initial_abundance"]]
  } else {
    cli_abort(
      c(
        "initial_abundance must be a numeric vector, raster, matrix, or
                array.",
        "x" = "initial_abundance is
        {.obj_type_friendly {inputs[['initial_abundance']]}}."
      )
    )
  }

  segments <- inputs[["segments"]] <-
    inputs[["stages"]]*inputs[["compartments"]]

  if (nrow(inputs[["initial_abundance"]]) != segments) {
    cli_abort(c(
      "initial_abundance has {nrow(inputs[['initial_abundance']])} row{?s}.",
      "x" = "There should be {segments} row{?s} (compartments x stages)."
    ))
  }
  if (ncol(inputs[["initial_abundance"]]) != inputs[["populations"]]) {
    cli_abort(c(
      "initial_abundance has {ncol(inputs[['initial_abundance']])} column{?s}.",
      "x" = "There should be {inputs[['populations']]} column{?s}."
    ))
  }

  # Check carrying capacity
  if (!is.null(inputs[["region"]]) &&
      inputs[["region"]][["use_raster"]] &&
      inherits(class(inputs[["carrying_capacity"]]),
               c("RasterLayer", "RasterStack", "RasterBrick"))) {
    inputs[["carrying_capacity_matrix"]] <- matrix(
      inputs[["carrying_capacity"]][inputs[["region"]][["region_indices"]]],
      nrow = inputs[["populations"]]
    )
  } else {
    inputs[["carrying_capacity_matrix"]] <- matrix(
      inputs[["carrying_capacity"]], nrow = inputs[["populations"]]
    )
  }
  if (anyNA(inputs[["carrying_capacity_matrix"]])) {
    inputs[["carrying_capacity_matrix"]][which(
      is.na(inputs[["carrying_capacity_matrix"]])
    )] <- 0
  }
  if (nrow(inputs[["carrying_capacity_matrix"]]) != inputs[["populations"]]) {
    cli_abort(c(
      "carrying_capacity has {nrow(inputs[['carrying_capacity_matrix']])} row{?s}.",
      "x" = "There should be {inputs[['populations']]} row{?s}."
    ))
  }
  if (!(ncol(inputs[["carrying_capacity_matrix"]]) %in%
        c(1, inputs[["time_steps"]]))) {
    cli_abort(c(
      "carrying_capacity has {ncol(inputs[['carrying_capacity_matrix']]}
      column{?s}.",
      "x" = "There should be 1 or {inputs[['time_steps']]} column{?s}."
    ))
  }
  if (any(inputs[["carrying_capacity_matrix"]] != floor(inputs[["carrying_capacity_matrix"]]))) {
    cli_abort(c(
      "carrying_capacity has values that are not whole numbers.",
      "i" = "You should ensure that carrying_capacity contains only whole numbers."
    ))
  }
  carrying_capacity_t_max <- inputs[["carrying_capacity_t_max"]] <-
    ncol(inputs[["carrying_capacity_matrix"]])
  if (carrying_capacity_t_max == 1) {
    # no temporal trend in K
    inputs[["carrying_capacity"]] <- inputs[["carrying_capacity_matrix"]][, 1]
  }

  # Check breeding season length
  if (!is.null(inputs[["breeding_season_length"]])) {
    if (!is.null(inputs[["region"]]) &&
    inputs[["region"]][["use_raster"]] &&
    any(c("RasterLayer", "RasterStack", "RasterBrick") %in%
    class(inputs[["breeding_season_length"]]))) {
      # Convert the RasterBrick to a matrix
      inputs[["breeding_season_matrix"]] <- raster::values(inputs[["breeding_season_length"]])

      # Subset the matrix
      inputs[["breeding_season_matrix"]] <- inputs[["breeding_season_matrix"]][inputs[["region"]][["region_indices"]],]
    } else {
      inputs[["breeding_season_matrix"]] <- matrix(
        inputs[["breeding_season_length"]], nrow = inputs[["populations"]]
      )
    }
    if (anyNA(inputs[["breeding_season_matrix"]])) {
      cli_abort(c("There are {sum(is.na(inputs[['breeding_season_matrix']]))}
                  missing values in the breeding_season_length object."))
    }
    breeding_season_t_max <- inputs[["breeding_season_t_max"]] <-
      ncol(inputs[["breeding_season_matrix"]])
    if (inputs[['seasons']] != 2) {
      cli_abort(c("breeding_season_length input only works in the two-season
                case (there is no way to infer the length of two non-breeding
                seasons from the length of one.)",
                "You have specified {inputs[['seasons']]} season{?s}."))
    }
    if (any(inputs[["breeding_season_matrix"]] != floor(inputs[["breeding_season_matrix"]]))) {
      cli_abort(c(
        "breeding_season has values that are not whole numbers.",
        "i" = "You should ensure that breeding_season contains only whole numbers."
      ))
    }
  } else if (!is.null(inputs[["season_lengths"]])) {
    if (length(inputs[["season_lengths"]]) != inputs[["seasons"]]) {
      cli_abort(c(
        "The length of `season_lengths` must equal `seasons`.",
        "i" = "`seasons` = {inputs[['seasons']]}.",
        "x" = "`season_lengths` = {inputs[['season_lengths']]}."
      ))
    }
  } else {
    inputs[["season_lengths"]] <- rep(
      round(365/inputs[['seasons']]), inputs[['seasons']]
    )
  }
  if (!is.null(inputs[["season_lengths"]]) &&
      !is.null(inputs[["breeding_season_matrix"]])) {
    cli_abort("You cannot specify both `season_lengths` and
              `breeding_season_length` inputs. Choose one.")
  }

  demography <- c("mortality", "mortality_unit", "mortality_mask",
    "fecundity", "fecundity_unit", "fecundity_mask", "transmission",
    "transmission_unit", "transmission_mask", "recovery", "recovery_unit",
    "recovery_mask", "dispersal")

  if (!is.null(inputs[["verbose"]])) {
    if (!is.logical(inputs[["verbose"]]) | length(inputs[["verbose"]]) > 1) {
      cli_abort(c("`verbose` must TRUE or FALSE."))
    }
  } else {
    inputs[["verbose"]] <- TRUE
  }

  for (element in demography) {
    if (is.vector(inputs[[element]]) && !is.null(names(inputs[[element]]))) {
      inputs[[element]] <- inputs[[element]][order(names(inputs[[element]]))]
      if (inputs[["verbose"]]) {
        cli_inform("Named lists of demographic rates have been ordered
               alphabetically.")
      }
    }
    if (is.list(inputs[[element]])) {
      for (i in 1:length(inputs[[element]])) {
        if (is.vector(inputs[[element]][[i]]) && !is.null(names(inputs[[element]][[i]]))) {
          inputs[[element]][[i]] <- inputs[[element]][[i]][order(names(inputs[[element]][[i]]))]
          if (inputs[["verbose"]]) {
            cli_inform("Named lists of demographic rates have been ordered
               alphabetically.")
          }
        }
      }
    }
  }

  # Check fecundity and survival
  mortality <- inputs[["mortality"]]
  if (is.list(mortality)) {
    if (length(mortality) != inputs[['seasons']]) {
      seasons <- inputs[['seasons']]
      cli_abort(c(
        "If {.var mortality} is given as a list, the list must be the same
        length as {.var seasons}.",
        "*" = "{.var mortality} is length {length(mortality)}.",
        "*" = "{.var seasons} is {seasons}."
      ))
    }
    if (any(map(mortality, length) != segments)) {
      cli_abort(
        c("Each vector inside the {.var mortality} list must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage & compartment.")
      )
    }
    if (mortality |> map_lgl(is.list) |> any()) {
      mortality <- inputs[["mortality"]] <- mortality |>
        map_if(is.list, flatten_dbl)
    }
  } else if (is.vector(mortality)) {
    if (length(mortality) != segments) {
      cli_abort(
        c("The {.var mortality} vector must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment.")
      )
    }
  } else {
    cli_abort(
      c("{.var mortality} must be a vector or list, not
        {.obj_type_friendly {mortality}}.")
    )
  }
  if (is.vector(mortality) && inputs[['seasons']] > 1 && !is.list(mortality)) {
    mortality <- inputs[["mortality"]] <- replicate(
      inputs[['seasons']], mortality, simplify = FALSE
    )
  }
  if (is.null(inputs[["mortality_unit"]])) {
    if (is.list(mortality)) {
      inputs[["mortality_unit"]] <- replicate(
        length(mortality),
        rep(0, length(mortality[[1]])),
        simplify = FALSE
      )
    } else {
      inputs[["mortality_unit"]] <- rep(0, segments)
    }
  } else {
    if (is.list(mortality) && !is.list(inputs[["mortality_unit"]])) {
      inputs[["mortality_unit"]] <- replicate(length(mortality),
                                              inputs[["mortality_unit"]],
                                              simplify = FALSE)
    }
    unit_values <- inputs[["mortality_unit"]] |> flatten_dbl() |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("mortality_unit values must be 0 or 1"))
    }
    if (length(unit_values)==1) {
      inputs[["mortality_unit"]] <- replicate(length(mortality),
                                              rep(unit_values, 8),
                                              simplify = FALSE)
    }
    if (length(mortality) != length(inputs[["mortality_unit"]])) {
      mortality_unit <- inputs[["mortality_unit"]]
      cli_abort(c("`mortality` and `mortality_unit` must be the same length.",
                  "*" = "`mortality` is length {length(mortality)}.",
                  "*" = "`mortality_unit` is length {length(mortality_unit)}."))
    }
    if (inputs[["mortality_unit"]] |> map_lgl(is.list) |> any()) {
      mortality_unit <- inputs[["mortality_unit"]] <- mortality_unit |>
        map_if(is.list, flatten_dbl)
    }
  }

  fecundity <- inputs[["fecundity"]]
  if (is.list(fecundity)) {
    if (length(fecundity) != inputs[['seasons']]) {
      seasons <- inputs[['seasons']]
      cli_abort(c(
        "If {.var fecundity} is given as a list, the list must be the same
        length as {.var seasons}.",
        "*" = "{.var fecundity} is length {length(fecundity)}.",
        "*" = "{.var seasons} is {seasons}."
      ))
    }
    if (any(lengths(fecundity) != segments)) {
      odd_indices <- which(lengths(fecundity) != segments)
      if (is.vector(inputs[["fecundity_mask"]])) {
        fecundity_mask <- inputs[["fecundity_mask"]]
        suppressWarnings(
          fecundity <- inputs[["fecundity"]] <- map2(fecundity[odd_indices],
                                                     fecundity_mask[odd_indices],
                                                     \(x, y) {
                                                       z <- rep(0, segments);
                                                       z[as.logical(y)] <- x})
        )
      }
    }
    if (any(map(fecundity, length) != segments)) {
      cli_abort(
        c("Each vector inside the {.var fecundity} list must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment. Or, provide a fecundity mask so
          that fecundity values may be assigned to the appropriate stages and
          compartments.")
      )
    }
    if (fecundity |> map_lgl(is.list) |> any()) {
      fecundity <- inputs[["fecundity"]] <- fecundity |>
        map_if(is.list, flatten_dbl)
    }
  } else if (is.vector(fecundity)) {
    if (length(fecundity) != segments) {
      if (is.vector(inputs[["fecundity_mask"]])) {
        fecundity_mask <- inputs[["fecundity_mask"]]
        z <- rep(0, segments)
        multiplier <- sum(fecundity_mask)/length(fecundity)
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
      c("{.var fecundity} must be a vector or list, not
        {.obj_type_friendly {fecundity}}.")
    )
  }
  if (is.vector(fecundity) && inputs[['seasons']] > 1 && !is.list(fecundity)) {
    fecundity <- inputs[["fecundity"]] <- replicate(
      inputs[['seasons']], fecundity, simplify = FALSE
    )
  }
  if (is.null(inputs[["fecundity_unit"]])) {
    if (is.list(fecundity)) {
      inputs[["fecundity_unit"]] <- replicate(
        length(fecundity),
        rep(0, length(fecundity[[1]])),
        simplify = FALSE
      )
    } else {
      inputs[["fecundity_unit"]] <- rep(0, segments)
    }
  } else {
    if (is.list(fecundity) && !is.list(inputs[["fecundity_unit"]])) {
      inputs[["fecundity_unit"]] <- replicate(length(fecundity),
                                              inputs[["fecundity_unit"]],
                                              simplify = FALSE)
    }
    unit_values <- inputs[["fecundity_unit"]] |> flatten_dbl() |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("fecundity_unit values must be 0 or 1. Unit values are
                  {unit_values}."))
    }
    if (length(fecundity) != length(inputs[["fecundity_unit"]])) {
      fecundity_unit <- inputs[["fecundity_unit"]]
      cli_abort(c("`fecundity` and `fecundity_unit` must be the same length.",
                  "*" = "`fecundity` is length {length(fecundity)}.",
                  "*" = "`fecundity_unit` is length {length(fecundity_unit)}."))
    }
    if (any(lengths(inputs[["fecundity_unit"]])==1)) {
      single_indices <- which(lengths(inputs[["fecundity_unit"]])==1)
      inputs[["fecundity_unit"]][single_indices] <- map(single_indices, \(x) {
        rep(inputs[["fecundity_unit"]][[x]], segments)
      })
    }
    if (any(lengths(fecundity) != lengths(inputs[["fecundity_unit"]]))) {
      fecundity_unit <- inputs[["fecundity_unit"]]
      cli_abort(c("vectors inside `fecundity` and `fecundity_unit` must be the
                  same length.",
                  "*" = "`fecundity` vectors are lengths {lengths(fecundity)}.",
                  "*" = "`fecundity_unit` vectors are lengths
                  {lengths(fecundity_unit)}."))
    }
    if (inputs[["fecundity_unit"]] |> map_lgl(is.list) |> any()) {
      inputs[["fecundity_unit"]] <- inputs[["fecundity_unit"]] |>
        map_if(is.list, list_c)
    }
  }
  if (is.null(inputs[["fecundity_mask"]])) {
    if(is.list(fecundity)) {
      inputs[["fecundity_mask"]] <- replicate(
        length(fecundity),
        rep(1, length(fecundity[[1]])),
        simplify = FALSE
      )
    } else {
      inputs[["fecundity_mask"]] <- rep(1, segments)
    }
  } else {
    if (is.list(fecundity) && !is.list(inputs[["fecundity_mask"]])) {
      if (length(inputs[["fecundity_mask"]]) == length(fecundity[[1]])) {
        inputs[["fecundity_mask"]] <- replicate(length(fecundity),
                                                inputs[["fecundity_mask"]],
                                                simplify = FALSE)
      } else {
        cli_abort(c("Vectors inside `fecundity` and `fecundity_mask` must be the
                    same length.",
                    "*" = "`fecundity` has vectors of length
                          {length(fecundity[[1]])}.",
                    "*" = "`fecundity_mask` has vectors of length
                          {length(inputs[['fecundity_mask']])}."))
      }
    }
    unit_values <- inputs[["fecundity_mask"]] |> flatten_dbl() |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("fecundity_mask values must be 0 or 1"))
    }
    if (is.list(fecundity) &&
        length(fecundity) != length(inputs[["fecundity_mask"]])) {
      fecundity_mask <- inputs[["fecundity_mask"]]
      cli_abort(c("`fecundity` and `fecundity_mask` must be the same length.",
                  "*" = "`fecundity` is length {length(fecundity)}.",
                  "*" = "`fecundity_mask` is length {length(fecundity_mask)}."))
    }
    if (inputs[["fecundity_mask"]] |> map_lgl(is.list) |> any()) {
      inputs[["fecundity_mask"]] <- inputs[["fecundity_mask"]] |>
        map_if(is.list, flatten_dbl)
    }
  }
  fecundity <- inputs[["fecundity"]] <- map2(fecundity,
                                             inputs[["fecundity_mask"]],
                     ~ ifelse(.y == 0, 0, .x))

  # Check transmission and recovery
  transmission <- inputs[["transmission"]]
  if (is.list(transmission)) {
    if (length(transmission) != inputs[['seasons']]) {
      seasons <- inputs[['seasons']]
      cli_abort(c(
        "If {.var transmission} is given as a list, the list must be the same
        length as {.var seasons}.",
        "*" = "{.var transmission} is length {length(transmission)}.",
        "*" = "{.var seasons} is {seasons}."
      ))
    }
    if (any(lengths(transmission) != segments)) {
      odd_indices <- which(lengths(transmission) != segments)
      if (is.list(inputs[["transmission_mask"]])) {
        transmission_mask <- inputs[["transmission_mask"]]
        suppressWarnings(
          transmission <- inputs[["transmission"]] <- map2(transmission[odd_indices],
                                                     transmission_mask[odd_indices],
                                                     \(x, y) {
                                                       z <- rep(0, segments);
                                                       z[as.logical(y)] <- x})
        )
      } else if (is.vector(inputs[["transmission_mask"]]) &&
                 length(inputs[["transmission_mask"]]) == segments) {
          transmission_mask <- inputs[["transmission_mask"]]
          transmission <- inputs[["transmission"]] <- map(
            transmission[odd_indices], \(x) {
              z <- rep(0, segments)
              z[as.logical(transmission_mask)] <- x
              return(z)
            }
          )
      }
    }
    if (any(map(transmission, length) != segments)) {
      cli_abort(
        c("Each vector inside the {.var transmission} list must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment. Or, provide a transmission mask so
          that transmission values may be assigned to the appropriate stages and
          compartments.")
      )
    }
    if (transmission |> map_lgl(is.list) |> any()) {
      transmission <- inputs[["transmission"]] <- transmission |>
        map_if(is.list, list_c)
    }
  } else if (is.vector(transmission)) {
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
  if (is.vector(transmission) && inputs[['seasons']] > 1 && !is.list(transmission)) {
    transmission <- inputs[["transmission"]] <- replicate(
      inputs[['seasons']], transmission, simplify = FALSE
    )
  }
  if (is.null(inputs[["transmission_unit"]])) {
    if (is.list(transmission)) {
      inputs[["transmission_unit"]] <- replicate(
        length(transmission),
        rep(0, length(transmission[[1]])),
        simplify = FALSE
      )
    } else {
      inputs[["transmission_unit"]] <- rep(0, segments)
    }
  } else {
    if (is.list(transmission) && !is.list(inputs[["transmission_unit"]])) {
      inputs[["transmission_unit"]] <- replicate(length(transmission),
                                                 inputs[["transmission_unit"]],
                                                 simplify = FALSE)
    }
    unit_values <- inputs[["transmission_unit"]] |> flatten_dbl() |> unique()
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
    if (any(lengths(inputs[["transmission_unit"]])==1)) {
      single_indices <- which(lengths(inputs[["transmission_unit"]])==1)
      inputs[["transmission_unit"]][single_indices] <- map(single_indices, \(x) {
        rep(inputs[["transmission_unit"]][[x]], segments)
      })
    }
    if (inputs[["transmission_unit"]] |> map_lgl(is.list) |> any()) {
      inputs[["transmission_unit"]] <- inputs[["transmission_unit"]] |>
      map_if(is.list, list_c)
    }
  }
  if (is.null(inputs[["transmission_mask"]])) {
    if(is.list(transmission)) {
      inputs[["transmission_mask"]] <- replicate(
        length(transmission),
        rep(0, segments) |> replace(1:stages, 1),
        simplify = FALSE
      )
    } else {
      inputs[["transmission_mask"]] <- rep(0, segments) |> replace(1:stages, 1)
    }
  } else {
    if (is.list(transmission) && !is.list(inputs[["transmission_mask"]])) {
      inputs[["transmission_mask"]] <- replicate(length(transmission),
                                                 inputs[["transmission_mask"]],
                                                 simplify = FALSE)
    }
    unit_values <- inputs[["transmission_mask"]] |> flatten_dbl() |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("transmission_mask values must be 0 or 1"))
    }
    if (is.list(transmission) &&
        length(transmission) != length(inputs[["transmission_mask"]])) {
      transmission_mask <- inputs[["transmission_mask"]]
      cli_abort(c("`transmission` and `transmission_mask` must be the same length.",
                  "*" = "`transmission` is length {length(transmission)}.",
                  "*" = "`transmission_mask` is length {length(transmission_mask)}."))
    }
    if (inputs[["transmission_mask"]] |> map_lgl(is.list) |> any()) {
      inputs[["transmission_mask"]] <- inputs[["transmission_mask"]] |>
        map_if(is.list, list_c)
    }
  }
  if (is.list(transmission) && length(transmission[[1]]) != segments &&
      is.list(inputs[['transmission_mask']])) {
    transmission <- inputs[["transmission"]] <- map2(transmission,
    inputs[["transmission_mask"]],
    function(x, mask) {
      z <- rep(0, length(mask))
      z[as.logical(mask)] <- x
      z
    })
  }
  apply_mask <- function(transmission_unit, transmission_mask) {
    if (is.list(transmission_unit) && is.list(transmission_mask)) {
      mapply(function(unit, mask) {
        unit_index <- 1
        result <- integer(length(mask))
        for (i in seq_along(mask)) {
          if (mask[i] == 1) {
            result[i] <- unit[unit_index]
            unit_index <- unit_index + 1
          } else {
            result[i] <- 0
          }
        }
        result
      }, transmission_unit, transmission_mask, SIMPLIFY = FALSE)
    } else if (is.vector(transmission_unit) && is.vector(transmission_mask)) {
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
    } else if (is.vector(transmission_unit) && is.list(transmission_mask)) {
      lapply(transmission_mask, function(mask) {
        unit_index <- 1
        result <- integer(length(mask))
        for (i in seq_along(mask)) {
          if (mask[i] == 1) {
            result[i] <- transmission_unit[unit_index]
            unit_index <- unit_index + 1
          } else {
            result[i] <- 0
          }
        }
        result
      })
    } else {
      cli_abort(c("Incongruity between mask vector and unit vector. Please
                  make sure that the sum of the mask vector matches the length
                  of the unit vector."))
    }
  }
  # use the function to apply the mask to the transmission unit
  inputs[["transmission_unit"]] <- apply_mask(inputs[["transmission_unit"]],
                                              inputs[["transmission_mask"]])
  if (any(lengths(transmission) != lengths(inputs[["transmission_unit"]]))) {
    transmission_unit <- inputs[["transmission_unit"]]
    cli_abort(c("vectors inside `transmission` and `transmission_unit` must be the
                  same length.",
                "*" = "`transmission` vectors are lengths {lengths(transmission)}.",
                "*" = "`transmission_unit` vectors are lengths
                  {lengths(transmission_unit)}."))
  }
  inputs[["transmission"]] <- map2(inputs[['transmission']],
                                 inputs[["transmission_mask"]],
                                 ~ ifelse(.y == 0, 0, .x))

  recovery <- inputs[["recovery"]]
  if (is.null(recovery)) {
    recovery <- inputs[["recovery"]] <- rep(0, segments)
  }
  if (is.list(recovery)) {
    if (length(recovery) != inputs[['seasons']]) {
      seasons <- inputs[['seasons']]
      cli_abort(c(
        "If {.var recovery} is given as a list, the list must be the same
        length as {.var seasons}.",
        "*" = "{.var recovery} is length {length(recovery)}.",
        "*" = "{.var seasons} is {seasons}."
      ))
    }
    if (any(lengths(recovery) != segments)) {
      odd_indices <- which(lengths(recovery) != segments)
      if (is.list(inputs[["recovery_mask"]])) {
        recovery_mask <- inputs[["recovery_mask"]]
        suppressWarnings(
          recovery <- inputs[["recovery"]] <- map2(recovery[odd_indices],
                                                     recovery_mask[odd_indices],
                                                     \(x, y) {
                                                       z <- rep(0, segments);
                                                       z[as.logical(y)] <- x})
        )
      } else if (is.vector(inputs[["recovery_mask"]]) &&
                 length(inputs[["recovery_mask"]]) == segments) {
          recovery_mask <- inputs[["recovery_mask"]]
          recovery <- inputs[["recovery"]] <- map(
            recovery[odd_indices], \(x) {
              z <- rep(0, segments)
              z[as.logical(recovery_mask)] <- x
              return(z)
            }
          )
      }
    }
    if (any(map(recovery, length) != segments)) {
      cli_abort(
        c("Each vector inside the {.var recovery} list must have
          {inputs[['stages']]*inputs[['compartments']]} elements, one for each
          combination of stage and compartment. Or, provide a recovery mask so
          that recovery values may be assigned to the appropriate stages and
          compartments.")
      )
    }
    if (recovery |> map_lgl(is.list) |> any()) {
      recovery <- inputs[["recovery"]] <- recovery |>
        map_if(is.list, list_c)
    }

  } else if (is.vector(recovery)) {
    if (length(recovery) != segments) {
      if (is.vector(inputs[["recovery_mask"]])) {
        recovery_mask <- inputs[["recovery_mask"]]
        z <- rep(0, segments)
        multiplier <- sum(recovery_mask)/length(recovery)
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
      c("{.var recovery} must be a vector or list, not
        {.obj_type_friendly {recovery}}.")
    )
  }
  if (is.vector(recovery) && inputs[['seasons']] > 1 && !is.list(recovery)) {
    recovery <- inputs[["recovery"]] <- replicate(
      inputs[['seasons']], recovery, simplify = FALSE
    )
  }
  if (is.null(inputs[["recovery_unit"]])) {
    if (is.list(recovery)) {
      inputs[["recovery_unit"]] <- replicate(
        length(recovery),
        rep(0, length(recovery[[1]])),
        simplify = FALSE
      )
    } else {
      inputs[["recovery_unit"]] <- rep(0, segments)
    }
  } else {
    if (is.list(recovery) && !is.list(inputs[["recovery_unit"]])) {
      inputs[["recovery_unit"]] <- replicate(length(recovery),
                                             inputs[["recovery_unit"]],
                                             simplify = FALSE)
    }
    unit_values <- inputs[["recovery_unit"]] |> flatten_dbl() |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("recovery_unit values must be 0 or 1. Unit values are
                  {unit_values}."))
    }
    if (length(recovery) != length(inputs[["recovery_unit"]])) {
      recovery_unit <- inputs[["recovery_unit"]]
      cli_abort(c("`recovery` and `recovery_unit` must be the same length.",
                  "*" = "`recovery` is length {length(recovery)}.",
                  "*" = "`recovery_unit` is length {length(recovery_unit)}."))
    }
    if (any(lengths(inputs[["recovery_unit"]])==1)) {
      single_indices <- which(lengths(inputs[["recovery_unit"]])==1)
      inputs[["recovery_unit"]][single_indices] <- map(single_indices, \(x) {
        rep(inputs[["recovery_unit"]][[x]], segments)
      })
    }
    if (inputs[["recovery_unit"]] |> map_lgl(is.list) |> any()) {
      inputs[["recovery_unit"]] <- inputs[["recovery_unit"]] |>
        map_if(is.list, list_c)
    }
  }
  if (is.null(inputs[["recovery_mask"]])) {
    if (inputs[["compartments"]] > 1) {
      recovery_mask <- rep(0, segments) |> replace((stages+1):(stages*2), 1)
    } else if (inputs[["compartments"]] == 1) {
      recovery_mask <- rep(0, segments)
    }
    if(is.list(recovery)) {
      inputs[["recovery_mask"]] <- replicate(
        length(recovery),
        recovery_mask,
        simplify = FALSE
      )
    } else {
      inputs[["recovery_mask"]] <- recovery_mask
    }
  } else {
    if (is.list(recovery) && !is.list(inputs[["recovery_mask"]])) {
      inputs[["recovery_mask"]] <- replicate(length(recovery),
                                             inputs[["recovery_mask"]],
                                             simplify = FALSE)
    }
    unit_values <- inputs[["recovery_mask"]] |> flatten_dbl() |> unique()
    if (!all(unit_values %in% c(0, 1))) {
      cli_abort(c("recovery_mask values must be 0 or 1"))
    }
    if (is.list(recovery) &&
        length(recovery) != length(inputs[["recovery_mask"]])) {
      recovery_mask <- inputs[["recovery_mask"]]
      cli_abort(c("`recovery` and `recovery_mask` must be the same length.",
                  "*" = "`recovery` is length {length(recovery)}.",
                  "*" = "`recovery_mask` is length {length(recovery_mask)}."))
    }
    if (inputs[["recovery_mask"]] |> map_lgl(is.list) |> any()) {
      inputs[["recovery_mask"]] <- inputs[["recovery_mask"]] |>
        map_if(is.list, list_c)
    }
  }
  # use the function to apply the mask to the recovery unit
  inputs[["recovery_unit"]] <- apply_mask(inputs[["recovery_unit"]],
                                          inputs[["recovery_mask"]])
  if (any(lengths(recovery) != lengths(inputs[["recovery_unit"]]))) {
    recovery_unit <- inputs[["recovery_unit"]]
    cli_abort(c("vectors inside `recovery` and `recovery_unit` must be the
                  same length.",
                "*" = "`recovery` vectors are lengths {lengths(recovery)}.",
                "*" = "`recovery_unit` vectors are lengths
                  {lengths(recovery_unit)}."))
  }
  recovery <- inputs[["recovery"]] <- map2(recovery,
                                           inputs[["recovery_mask"]],
                    ~ ifelse(.y == 0, 0, .x))


  if (is.null(inputs[["density_stages"]])) { # default is all
    inputs[["density_stages"]] <- rep(1, inputs[["stages"]])
  }

  # Dispersal parameters
  dispersal_stages <- inputs[["dispersal_stages"]]
  if (is.null(dispersal_stages)) {
    inputs[["dispersal_stages"]] <- array(1, inputs[["stages"]])
  }
  dispersal_source_n_k <- inputs[["dispersal_source_n_k"]]
  if (is.null(dispersal_source_n_k)) {
    # try aliases
    inputs[["dispersal_source_n_k"]] <- list(
      cutoff = inputs[["dispersal_source_n_k_cutoff"]],
      threshold = inputs[["dispersal_source_n_k_threshold"]]
    )
  }
  dispersal_target_k <- inputs[["dispersal_target_k"]]
  if (is.null(dispersal_target_k)) { # try alias
    inputs[["dispersal_target_k"]] <- inputs[["dispersal_k_threshold"]]
  }
  dispersal_target_n <- inputs[["dispersal_target_n"]]
  if (is.null(dispersal_target_n)) {
    # try aliases
    inputs[["dispersal_target_n"]] <- list(
      threshold = inputs[["dispersal_n_threshold"]],
      cutoff = inputs[["dispersal_n_cutoff"]]
    )
  }
  dispersal_target_n_k <- inputs[["dispersal_target_n_k"]]
  if (is.null(dispersal_target_n_k)) {
    # try aliases
    inputs[["dispersal_target_n_k"]] <- list(
      threshold = inputs[["dispersal_target_n_k_threshold"]],
      cutoff = inputs[["dispersal_target_n_k_cutoff"]]
    )
  }
  if (!is.null(inputs[["dispersal"]]) &&
      !length(inputs[["dispersal"]]) %in% c(1, stages, compartments, segments)) {
    cli_abort(c('`dispersal` must be a list of length 1 or
                {inputs[["stages"]]} (number of stages),
                {compartments} (number of compartments), or
                {segments} (stages * compartments).",
                "x" = "{.var dispersal} is length {length(inputs[["dispersal"]])}.'))
  }

  if(is.null(inputs[["dispersal_type"]]) && !is.null(inputs[["dispersal"]])) {
      inputs[["dispersal_type"]] <- "pooled"
  }

  if (!is.null(inputs[["dispersal_type"]]) && !is.null(inputs[["dispersal"]])) {
      if (inputs[["dispersal_type"]] == "pooled") {
          if (length(inputs[["dispersal"]]) != 1) {
              cli_abort(c("When `dispersal_type` is pooled, `dispersal` must be length 1."))
          }
      } else if (inputs[["dispersal_type"]] == "stages" && !is.null(inputs[["stages"]])) {
          if (length(inputs[["dispersal"]]) != inputs[["stages"]]) {
              cli_abort("Length of `dispersal` must be equal to the number of stages.")
          }
      } else if (inputs[["dispersal_type"]] == "compartments" && !is.null(inputs[["compartments"]])) {
          if (length(inputs[["dispersal"]]) != inputs[["compartments"]]) {
              cli_abort("Length of `dispersal` must be equal to the number of compartments.")
          }
      } else if (inputs[["dispersal_type"]] == "segments" && exists("segments")) {
          if (length(inputs[["dispersal"]]) != segments) {
              cli_abort("Length of `dispersal` must be equal to `stages*compartments`.")
          }
      } else {
          cli_abort("Invalid value for `dispersal_type`.")
      }
  }

  # Abundance threshold
  abundance_threshold <- inputs[["abundance_threshold"]]
  if (is.null(abundance_threshold)) {
    inputs[["abundance_threshold"]] <- abundance_threshold <- 0
  }
  if (length(abundance_threshold) > 1) {
    cli_abort(
      c("{.var abundance_threshold} must be a single value.",
        "x" = "{.var abundance_threshold} is length
        {length(abundance_threshold)}.")
    )
  }

  # Simulation order
  simulation_order <- inputs[["simulation_order"]]
  if (is.vector(simulation_order) && !is.list(simulation_order)) {
    simulation_order <- inputs[["simulation_order"]] <- replicate(
      inputs[["seasons"]], inputs[['simulation_order']], simplify = FALSE
    )
  }
  if (!all(simulation_order |> flatten() |> map_lgl(is.character))) {
    cli_abort(c("{.var simulation_order} must contain only strings."))
  }
  if (!"results" %in% flatten(simulation_order)) {
    no_result_index <- simulation_order |> map("results") |> map_lgl(is_null)
    inputs[['simulation_order']][no_result_index] <-
      inputs[['simulation_order']][no_result_index] |> map(append, "results")
  }

  # Results selection and breakdown
  if (is.null(inputs[["results_selection"]])) {
    inputs[['results_selection']] <- c("abundance", "summarize")
  }
  results_selection <- inputs[["results_selection"]]
  if (!length(intersect(results_selection, c("abundance", "ema", "extirpation",
                                             "extinction_location", "harvested",
                                             "occupancy", "summarize",
                                             "replicate")))) {
    cli_abort(
      c('{.var results_selection} must contain at least one of the following:
        {c("abundance", "ema", "extirpation", "extinction_location",
        "harvested", "occupancy", "summarize", "replicate")}')
    )
  }

  if (is.null(inputs[["results_breakdown"]])) {
    inputs[["results_breakdown"]] <- "segments"
  }
  results_breakdown <- inputs[["results_breakdown"]]
  if (length(results_breakdown) > 1) {
    cli_abort(
      c("{.var results_breakdown} must be a character vector of length 1.",
        "x" = "{.var results_breakdown} is length {length(results_breakdown}.")
    )
  }

  if (!results_breakdown %in% c("pooled", "stages", "compartments", "segments")) {
    cli_abort(
      c("{.var results_breakdown} must be 'pooled', 'stages', 'compartments' or
        'segments.'",
        "x" = "{.var results_breakdown} is {results_breakdown}.")
    )
  }

  harvest <- inputs[["harvest"]]
  if (!is.null(harvest)) {
    if (!is.list(harvest) && !is.function(harvest)) {
      cli_abort(
        c("{.var harvest} must be a function or a list.",
          "x" = "{.var harvest} is {class(harvest)}.")
      )
    }
  }

  translocation <- inputs[["translocation"]]
  if (!is.null(translocation)) {
    if (!is.list(translocation) && !is.function(translocation)) {
      cli_abort(
        c("{.var translocation} must be a function or a list.",
          "x" = "{.var translocation} is {class(translocation)}.")
      )
    }
  }

  mortality_function <- inputs[["mortality_function"]]
  if (!is.null(mortality_function)) {
    if (!is.list(mortality_function) && !is.function(mortality_function)) {
      cli_abort(
        c("{.var mortality_function} must be a function or a list.",
          "x" = "{.var mortality_function} is {class(mortality_function)}.")
      )
    }
  }

  return(inputs)
}
