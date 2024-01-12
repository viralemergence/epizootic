#'Stage-based seasonal spatially explicit population-level disease model.
#'
#'Simulates a stage-based demographic population model and returns simulation
#'results across multiple replicate runs. Processes run at each simulation
#'time-step include:
#' \enumerate{
#'   \item Density dependence calculations (ceiling, logistic, or user-defined)
#'   \item Stage transition (stochastic) calculations
#'   \item Population growth/decline calculations
#'   \item Disease outbreak according to a compartmental model
#'   \item Dispersal calculations (default or user-defined)
#'   \item Results collection
#' }
#'Note that the breeding season is
#'always treated as the first season.
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
#'  \item{\code{region}}{A \code{\link[poems:Region]{Region}} object
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
#'  as per \code{\link{SpatialCorrelation}} class attributes.}
#'  \item{\code{mortality}}{A vector of mortality rates, one for each
#'  combination of stages and compartments. Assumed by default to be daily
#'  mortality rates unless indicated otherwise (see below). If mortality varies
#'  by season, a list of mortality vectors with the same length as `seasons`
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
#'  \item{\code{abundance_threshold}}{A quasi-extinction threshold below which a
#'  population becomes extinct. Default: 0.}
#'  \item{\code{demographic_stochasticity}}{Boolean for choosing demographic
#'  stochasticity for transition, dispersal, and/or other processes (default is
#'  TRUE).}
#'  \item{\code{transmission}}{A vector of transmission rates, one for each
#'  combination of stages and compartment for which transmission applies (see
#'  \code{transmission_mask} below. If transmission varies by season, a list of
#'  transmission vectors with the same length as `seasons` may be provided
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
#'  compartment can recover.}
#'  \item{\code{density_dependence}}{Selects a density dependence function, if
#'  desired. Choices are "none" (the default), "ceiling", "logistic" (Ricker),
#'  or a user-defined function (optionally nested in a list with additional
#'  attributes) for adjusting transition rates: \code{function(params)},
#'  where \code{params} is a list passed to the function containing:
#'       \describe{
#'         \item{\code{transition_array}}{3D array of transition rates: stages
#'         by stages by populations.}
#'         \item{\code{fecundity_mask}}{Matrix of 0-1 to indicate which
#'         (proportions) of transition rates refer to fecundity.}
#'         \item{\code{fecundity_max}}{Maximum transition fecundity rate
#'         (in Leslie/Lefkovitch matrix).}
#'         \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'         each population.}
#'         \item{\code{segment_abundance}}{Matrix of abundances for each stage-
#'         compartment combination (rows) and population (columns).}
#'         \item{\code{population_abundance}}{Array of summed population
#'         abundances for all stages.}
#'         \item{\code{density_abundance}}{Array of summed population abundances
#'         for stages affected by density.}
#'         \item{\code{growth_rate_max}}{Maximum growth rate value or array for
#'         populations.}
#'         \item{\code{occupied_indices}}{Array of indices for populations
#'         occupied at (current) time step.}
#'         \item{\code{calculate_multipliers}}{Function
#'         (\code{function(growth_rates)}) for finding multipliers (when stages
#'         \> 1) to apply to affected transitions that result in target growth
#'          rates (dominant eigenvalues).}
#'         \item{\code{apply_multipliers}}{Function
#'         (\code{function(transition_array, multipliers}) for applying
#'         multipliers (when stages \> 1) to the affected transition rates
#'         within a transition array (returns multiplied array).}
#'         \item{\code{simulator}}{\code{\link{SimulatorReference}} object with
#'         dynamically accessible \emph{attached} and \emph{results} lists.}
#'         \item{\code{optional attributes}}{Additional numeric attributes when
#'         density dependence is optionally nested in a list.}
#'       }
#'       and returns a transformed transition 3D array.
#'     }
#'  \item{\code{growth_rate_max}}{Maximum growth rate (utilized by density
#'  dependence processes).}
#'  \item{\code{density_stages}}{Array of booleans or numeric
#'  (0,1) for each stage to indicate which stages are affected by density
#'  (default is all).}
#'  \item{\code{translocation}}{An optional user-defined function (optionally
#'     nested in a list with additional attributes) for applying translocation
#'     or spatio-temporal management (to abundances): \code{function(params)},
#'     where \emph{params} is a list passed to the function containing:
#'       \describe{
#'         \item{\code{replicates}}{Number of replicate simulation runs.}
#'         \item{\code{time_steps}}{Number of simulation time steps.}
#'         \item{\code{years_per_step}}{Number of years per time step.}
#'         \item{\code{populations}}{Number of populations.}
#'         \item{\code{stages}}{Number of lifecycle stages.}
#'         \item{\code{demographic_stochasticity}}{Boolean for optionally
#'         choosing demographic stochasticity for the transformation.}
#'         \item{\code{density_stages}}{Array of booleans or numeric (0,1) for
#'         each stage to indicate which stages are affected by density.}
#'         \item{\code{r}}{Simulation replicate.}
#'         \item{\code{tm}}{Simulation time step.}
#'         \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'         each population at time step.}
#'         \item{\code{segment_abundance}}{Matrix of current abundance for each
#'         stage-compartment combination (rows) and population (columns) at time
#'         step.}
#'         \item{\code{occupied_indices}}{Array of indices for populations
#'         occupied at (current) time step.}
#'         \item{\code{simulator}}{\code{\link{SimulatorReference}} object with
#'         dynamically accessible \emph{attached} and \emph{results} lists.}
#'         \item{\code{additional attributes}}{Additional attributes when the
#'         transformation is optionally nested in a list.}
#'       }
#'       and returns a transformed stage abundance matrix (or a list with stage
#'       abundance and carrying capacity).
#'     }
#'  \item{\code{harvest}}{An optional user-defined function (optionally nested
#'  in a list with additional attributes) for applying harvesting (to
#'  abundances): \code{function(params)} as per translocation.}
#'  \item{\code{mortality_function}}{An optional user-defined function (optionally nested
#'  in a list with additional attributes) for applying mortality (to
#'  abundances): \code{function(params)} as per translocation.}
#'  \item{\code{dispersal}}{A list that is either length 1 or the same length as
#'  \code{stages}. If it is length 1, the same dispersal will be applied across
#'  all stages. Within each element of the list, there should be either a matrix
#'  of dispersal rates between populations (source columns to target rows) or a
#'  list of data frames of non-zero dispersal rates and indices for constructing
#'  a compact dispersal matrix, and optional changing rates over time (as per
#'  class \code{\link{DispersalGenerator}} \emph{dispersal_data} attribute).}
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
#'          \item{\code{r}}{Simulation replicate.} \item{\code{tm}}{Simulation
#'          time step.}
#'          \item{\code{carrying_capacity}}{Array of carrying capacity values
#'          for each population at time step.}
#'          \item{\code{segment_abundance}}{Matrix of abundance for each
#'          combination of stage and compartment (rows) and population (columns)
#'           at time step.}
#'          \item{\code{occupied_indices}}{Array of indices for populations
#'          occupied at time step.}
#'          \item{\code{simulator}}{\code{\link{SimulatorReference}} object
#'          with dynamically accessible \emph{attached} and \emph{results}
#'          lists.}
#'          \item{\code{additional attributes}}{Additional attributes when the
#'          transformation is optionally nested in a list.}
#'          }
#'       and returns a transformed stage abundance matrix.
#'  }
#'  \item{\code{simulation_order}}{A list the same length as \code{seasons}.
#'  Each element in the list is a vector of named simulation processes in the
#'  desired order. Processes must be one of "transition", "translocation",
#'  "harvest", "mortality", "dispersal", "season_functions", or "results."
#'  "season_functions" will be matched to the appropriate season (i.e., if
#'  "season_functions" appears in element 1 of the list, `season_functions[[1]]`
#'  will be called.) If the simulation processes are the same across seasons,
#'  then a single character vector may be provided. Required input.}
#'  \item{\code{additional transformation functions}}{Additional user-defined
#'  abundance transformation functions (optionally nested in lists with
#'  additional attributes) are utilised when listed in \emph{simulation_order}
#'  (function as per translocation).}
#'  \item{\code{results_selection}}{List of results selection from: "abundance"
#'  (default), "ema", "extirpation", "extinction_location", "harvested",
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
#'@return Selected simulation results as a nested list summarized (mean, sd,
#'  min, max) across multiple replicates (default), or 2-3D arrays including
#'  results for each replicate:
#' \describe{
#'   \item{\code{abundance}}{Matrix or 3D array of simulation abundance:
#'   \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'   deep).}
#'   \item{\code{abundance_stages}}{List of matrices or 3D arrays of simulation
#'   abundance for unique stage-compartment combinations when present: each
#'   \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'   deep).}
#'   \item{\code{all$abundance}}{Array or matrix of total abundance across
#'   populations: \emph{time_steps} (rows by \emph{replicates} columns).}
#'   \item{\code{all$abundance_stages}}{List of arrays or matrices of total
#'   abundance across populations for unique stage-compartment combinations when
#'   present: each \emph{time_steps} (rows by \emph{replicates} columns).}
#'   \item{\code{all$ema}}{Array of expected minimum abundance at each time step
#'   (averaged across replicates).}
#'   \item{\code{extirpation}}{Array or matrix of extirpation times:
#'   \emph{populations} (rows by \emph{replicates} columns).}
#'   \item{\code{all$extirpation}}{Array of extirpation time across populations
#'   for each replicate.}
#'   \item{\code{all$extinction_location}}{The weighted centroid of cells
#'   occupied in the time-step prior to the extirpation of all populations
#'   (if it occurred) for each replicate.}
#'   \item{\code{harvested}}{Matrix or 3D array of individuals harvested:
#'   \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'   deep).}
#'   \item{\code{harvested_stages}}{List of matrices or 3D arrays of individuals
#'   harvested for unique stage-compartment combinations when present: each
#'   \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'   deep).}
#'   \item{\code{all$harvested}}{Array or matrix of individuals harvested across
#'   populations: \emph{time_steps} (rows by \emph{replicates} columns).}
#'   \item{\code{all$harvested_stages}}{List of arrays or matrices of
#'   individuals harvested across populations for unique stage-compartment
#'   combinations when present: each \emph{time_steps} (rows by
#'   \emph{replicates} columns).}
#'   \item{\code{all$occupancy}}{Array or matrix of the number of populations
#'   occupied at each time-step: \emph{time_steps} (rows by \emph{replicates}
#'   columns).}
#'   \item{\code{additional results}}{Additional results may be attached via
#'   user-defined functions (using \code{params$simulator$results}).}
#' }
#'
#' @import poems
#' @import cli
#' @include check_simulator_inputs.R
#' @export disease_simulator

disease_simulator <- function(inputs) {

  # Check that all inputs are valid
  inputs <- check_simulator_inputs(inputs)
  list2env(inputs, envir = environment())

  population_abundance <- colSums(initial_abundance)

  # Simulator reference object for dynamic attachments and results
  # (accessed via user-defined functions)
  simulator <- SimulatorReference$new()

  # Generate simulation functions
  transition_function <- disease_transitions(stages, compartments)

  if (density_dependence != "none") {
    density_function <- population_density(
      populations,
      stage_matrix,
      fecundity_mask,
      fecundity_max,
      density_dependence,
      inputs$growth_rate_max,
      inputs$density_affects,
      density_stages,
      density_precision = inputs$density_precision,
      simulator
    )
    if ((is.function(density_dependence) ||
         is.list(density_dependence))) {
      density_dependence <- "adjusts_transitions"
    }
  }

  translocation_function <- population_transformation(
    replicates = replicates,
    time_steps = time_steps,
    populations = populations,
    demographic_stochasticity = demographic_stochasticity,
    density_stages = density_stages,
    transformation = inputs[["translocation"]],
    simulator = simulator,
    name = "translocation"
  )
  harvest_function <- population_transformation(
    replicates = replicates,
    time_steps = time_steps,
    populations = populations,
    demographic_stochasticity = demographic_stochasticity,
    density_stages = density_stages,
    transformation = inputs[["harvest"]],
    simulator = simulator,
    name = "harvest"
  )

  if ("mortality" %in% simulation_order) {
    mortality_function <- population_transformation(
      replicates = replicates,
      time_steps = time_steps,
      populations = populations,
      demographic_stochasticity = demographic_stochasticity,
      density_stages = density_stages,
      transformation = inputs[["mortality_function"]],
      simulator = simulator,
      name = "mortality"
    )
  }

  if (length(inputs[["dispersal"]]) == 1) {
    dispersal_function <- population_dispersal(
      replicates,
      time_steps,
      populations,
      demographic_stochasticity,
      density_stages,
      flatten(inputs[["dispersal"]]),
      dispersal_stages,
      dispersal_source_n_k,
      dispersal_target_k,
      dispersal_target_n,
      dispersal_target_n_k,
      simulator
    )
  } else if (length(inputs[["dispersal"]]) > 1) {
    dispersal_functions = list()
    for (d in 1:length(inputs[["dispersal"]])) {
      dispersal_functions[[d]] <- population_dispersal(
        replicates,
        time_steps,
        populations,
        demographic_stochasticity,
        density_stages,
        inputs[["dispersal"]][[d]],
        dispersal_stages,
        dispersal_source_n_k,
        dispersal_target_k,
        dispersal_target_n,
        dispersal_target_n_k,
        simulator
      )
    }
  }

  if (exists("season_functions")) {
    season_function_list <- list()
    for (i in 1:length(season_functions)) {
      season_function_list[[i]] <- disease_transformation(
        replicates,
        time_steps,
        seasons,
        compartments,
        populations,
        demographic_stochasticity,
        stages,
        abundance_threshold,
        mortality[[i]],
        mortality_unit[[i]],
        fecundity[[i]],
        fecundity_unit[[i]],
        fecundity_mask[[i]],
        transmission[[i]],
        transmission_unit[[i]],
        transmission_mask[[i]],
        recovery[[i]],
        recovery_unit[[i]],
        recovery_mask[[i]],
        transformation = season_functions[[i]],
        simulator,
        name = paste0("season", i, collapse = "_")
      )
    }
  }

  result_functions <- disease_results(
    replicates = replicates,
    time_steps = time_steps,
    seasons = seasons,
    stages = stages,
    compartments = compartments,
    coordinates = inputs[["coordinates"]],
    initial_abundance = initial_abundance,
    results_selection = results_selection,
    results_breakdown = results_breakdown
  )
  results_list <- result_functions$initialize_attributes()

  ### Replicates ###
  for (r in 1:replicates) {

    # Initialize populations
    segment_abundance <- initial_abundance
    population_abundance <- .colSums(initial_abundance,
                                     m = segments, n = populations)
    occupied_indices <- which(as.logical(population_abundance))
    occupied_populations <- length(occupied_indices)

    # Initialize harvested
    if ("harvested" %in% results_selection) {
      harvested <- array(0 , c(segments, populations))
    } else {
      harvested <- NULL
    }

    # (Re-)Initialize result collection variables
    results_list <- result_functions$initialize_replicate(results_list)

    ### Simulation time steps ###
    for (tm in 1:time_steps) {

      # Load carrying capacity for each population for time if there is a
      # temporal trend in K
      if (carrying_capacity_t_max > 1) {
        carrying_capacity <-
          carrying_capacity_matrix[, min(tm, carrying_capacity_t_max)]
      }

      ## Run simulation processes in configured order ##
      for (season in 1:seasons) {
        sim_order <- simulation_order[[season]]

        for (process in sim_order) {

          if (process == "transition") {

            if (occupied_populations) {
              # Perform stage-based transitions
              segment_abundance <- transition_function(segment_abundance,
                                                       occupied_indices)

              # Limit abundances to carrying capacity when "ceiling" density
              # dependence
              if (density_dependence == "ceiling") {
                segment_abundance <- density_function(carrying_capacity,
                                                      segment_abundance)
              }
            }
          }

          ## Translocation calculations ##
          if (process == "translocation" &&
              is.function(translocation_function)) {
            transformed <- translocation_function(r, tm, carrying_capacity,
                                                  segment_abundance,
                                                  occupied_indices)
            segment_abundance <- transformed[["segment_abundance"]]
            if ("carrying_capacity" %in% names(transformed)) {
              carrying_capacity <- transformed[["carrying_capacity"]]
            }
          }

          ## Harvest calculations ##
          if (process == "harvest") {
            if (occupied_populations && is.function(harvest_function)) {
              preharvest_abundance <- segment_abundance
              transformed <- harvest_function(r, tm, carrying_capacity,
                                              segment_abundance,
                                              occupied_indices)
              segment_abundance <- transformed$segment_abundance
              if ("carrying_capacity" %in% names(transformed)) {
                carrying_capacity <- transformed$carrying_capacity
              }
              harvested <- preharvest_abundance - segment_abundance
            } else {
              harvested <- 0*segment_abundance
            }
            if ("harvested" %in% results_selection) {
              results_list <- result_functions$calculate_at_timestep(r, tm,
                                                                     NULL,
                                                                     harvested,
                                                                     results_list)
            }
          }

          ## Mortality calculations ##
          if (occupied_populations && process == "mortality" &&
              is.function(mortality_function)) {
            transformed <- mortality_function(r, tm, carrying_capacity,
                                              segment_abundance,
                                              occupied_indices)
            segment_abundance <- transformed$segment_abundance
            if ("carrying_capacity" %in% names(transformed)) {
              carrying_capacity <- transformed$carrying_capacity
            }
          }

          ## Dispersal calculations ##
          if (occupied_populations && process == "dispersal" &&
              !is.null(dispersal)) {
            if (!is.null(dispersal_function)) {
              segment_abundance <- dispersal_function(r, tm, carrying_capacity,
                                                      segment_abundance,
                                                      occupied_indices)
            } else if (!is.null(dispersal_functions)) {
              indices <- map(1:stages, \(s) seq(s, segments, stages))
              dispersal_function <- dispersal_functions[[s]]
              stage_dispersal <- map(1:stages, \(s) {
                abundance <- dispersal_function(r, tm, carrying_capacity,
                                                segment_abundance[indices[[s]],],
                                                occupied_indices)
                return(abundance)
              })
              for (i in 1:stages) {
                segment_abundance[indices[[i]],] <- stage_dispersal[[i]]
              }
            }
          }

          if (process == "season_functions" && is.list(season_function_list)) {
            transformed <- season_function_list[[season]](r, tm,
                                                          carrying_capacity,
                                                          segment_abundance,
                                                          breeding_season_length,
                                                          occupied_indices)
            segment_abundance <- transformed$segment_abundance
            if ("carrying_capacity" %in% names(transformed)) {
              carrying_capacity <- transformed$carrying_capacity
            }
          }

          if (process == "results") {
            results_list <- result_functions$calculate_at_season(r, tm, season,
                                                                 segment_abundance,
                                                                 NULL,
                                                                 results_list)
          }
        } # End simulation order loop
      } # End season loop
    } # End time step loop

    results_list <- result_functions$calculate_at_replicate(r,
                                                            segment_abundance,
                                                            results_list)

  } # End replicate loop

  ## Finalize results calculation and collation ##
  results_list <- result_functions$finalize_attributes(results_list)

  return(c(results_list, simulator$results))

}
