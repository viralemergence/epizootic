#'Stage-based seasonal spatially explicit SIRI model of the
#'*Mycoplasma gallisepticum* outbreak in *Haemorhous mexicanus*.
#'
#'Simulates a stage-based demographic population model and returns simulation
#'results across multiple replicate runs. Processes run at each simulation
#'time-step include:
#' \enumerate{
#'   \item Density dependence calculations (ceiling, logistic, or user-defined)
#'   \item Environmental stochasticity calculations
#'   \item Stage transition (stochastic) calculations
#'   \item Population growth/decline calculations
#'   \item Disease outbreak according to a SIRI model
#'   \item Dispersal calculations (default or user-defined)
#'   \item Results collection
#' }
#'This simulator is not ready for use in any system other than the *Mycoplasma
#'gallisepticum* outbreak in the house finch. Note that the breeding season is
#'always treated as the first season.
#'
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
#'  If there is only one stage and one compartment and a region object is
#'  attached, then initial abundance may be provided in the form of a raster
#'  with the same specs as the region raster. Similarly, if there is only one
#'  stage/compartment combination you may provide a vector with length
#'  \code{populations}. Required input.}
#'  \item{\code{carrying_capacity}}{Array (matrix) of carrying capacity values
#'  at each population cell (\code{populations} rows by \code{time_steps}
#'  columns when across time). Required input.}
#'  \item{\code{breeding_season_length}}{Array (matrix) of breeding season
#'  length values in days at each population cell (\code{populations} rows by
#'  \code{time_steps} columns when across time).}
#'  \item{\code{season_lengths}}{Vector of season lengths in days. Length must
#'  equal \code{seasons}.}
#'  \item{\code{correlation}}{List containing either an environmental
#'  correlation matrix (correlation_matrix), a pre-calculated transposed
#'  (Cholesky) decomposition matrix (t_decomposition_matrix), or a compact
#'  transposed (Cholesky) decomposition matrix (t_decomposition_compact_matrix)
#'  and a corresponding map of population indices (t_decomposition_compact_map),
#'  as per \code{\link{SpatialCorrelation}} class attributes.}
#'  \item{\code{mortality}}{A vector of mortality rates, one for each combination
#'  of stages and compartments. Assumed by default to be daily mortality rates.
#'  Required input.}
#'  \item{\code{mortality_unit}}{A vector indicating whether mortality rates are
#'  daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all 0.}
#'  \item{\code{fecundity}}{A vector of fecundity rates, one for each
#'  combination of stages and compartments for which fecundity applies (see
#'  \code{fecundity_mask} below). Required input.}
#'  \item{\code{fecundity_unit}}{A vector indicating whether mortality rates are
#'  daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all 0.}
#'  \item{\code{fecundity_mask}}{A vector indicating which stages and
#'  compartments reproduce. Must be the same length as
#'  \code{stages * compartments}.}
#'  \item{\code{abundance_threshold}}{A quasi-extinction threshold below which a
#'  population becomes extinct.}
#'  \item{\code{demographic_stochasticity}}{Boolean for choosing demographic
#'  stochasticity for transition, dispersal, and/or other processes (default is
#'  TRUE).}
#'  \item(\code{transmission}}{A vector of transmission rates, one for each
#'  combination of stages and compartment for which transmission applies (see
#'  \code{transmission_mask} below. Required input.}
#'  \item{\code{transmission_unit}}{A vector indicating whether mortality rates
#'  are daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default: all
#'  0.}}
#'  \item{\code{transmission_mask}}{A vector indicating which compartments are
#'  subject to transmission (i.e., classes susceptible to infection.) Must be
#'  the same length as \code{compartments}.}
#'  \item{\code{recovery}}{A vector of recovery rates, one for each
#'  combination of stages and compartment for which recovery applies (see
#'  \code{recovery_mask} below.)}}
#'  \item{\code{recovery_unit}}{A vector
#'  indicating whether mortality rates are daily or seasonal. 1 indicates
#'  seasonal, 0 indicates daily. Default: all 0.}}
#'  \item{\code{recovery_mask}}{A
#'  vector indicating which compartments are subject to recovery (i.e., infected
#'  classes that can recover.) Must be the same length as \code{compartments}.}
#'  \item{\code{density_dependence}}{Selects a density dependence function, if
#'  desired. Choices are "ceiling", "logistic" (Ricker), or a user-defined
#'  function (optionally nested in a list with additional attributes) for
#'  adjusting transition rates: \code{function(params)}, where \code{params} is
#'  a list passed to the function containing:
#'       \describe{
#'         \item{\code{transition_array}}{3D array of transition rates: stages
#'         by stages by populations.}
#'         \item{\code{fecundity_mask}}{Matrix of 0-1 to indicate which
#'         (proportions) of transition rates refer to fecundity.}
#'         \item{\code{fecundity_max}}{Maximum transition fecundity rate
#'         (in Leslie/Lefkovitch matrix).}
#'         \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'         each population.}
#'         \item{\code{stage_abundance}}{Matrix of abundances for each stage-
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
#'         > 1) to apply to affected transitions that result in target growth
#'          rates (dominant eigenvalues).}
#'         \item{\code{apply_multipliers}}{Function
#'         (\code{function(transition_array, multipliers}) for applying
#'         multipliers (when stages > 1) to the affected transition rates within
#'         a transition array (returns multiplied array).}
#'         \item{\code{simulator}}{\code{\link{SimulatorReference}} object with
#'         dynamically accessible \emph{attached} and \emph{results} lists.}
#'         \item{\code{optional attributes}}{Additional numeric attributes when
#'         density dependence is optionally nested in a list.}
#'       }
#'       and returns a transformed transition 3D array.
#'     }
#'  \item{\code{growth_rate_max}}{Maximum growth rate (utilized by density
#'  dependence processes).}
#'  \item{\code{density_affects}}{Matrix of booleans or
#'  numeric (0-1) indicating the transition vital rates affected by density
#'  (default is all).}
#'  \item{\code{density_stages}}{Array of booleans or numeric
#'  (0,1) for each stage to indicate which stages are affected by density
#'  (default is all).}
#'  \item{\code{density_precision}}{Numeric precision of the
#'  calculated multipliers (used when stages > 1) applied to affected transition
#'  rates (default is 3 decimal places).}
#'  \item{\code{translocation}}{An optional user-defined function (optionally
#'     nested in a list with additional attributes) for applying translocation or
#'     spatio-temporal management (to abundances): \code{function(params)}, where
#'     \emph{params} is a list passed to the function containing:
#'       \describe{
#'         \item{\code{replicates}}{Number of replicate simulation runs.}
#'         \item{\code{time_steps}}{Number of simulation time steps.}
#'         \item{\code{years_per_step}}{Number of years per time step.}
#'         \item{\code{populations}}{Number of populations.}
#'         \item{\code{stages}}{Number of lifecycle stages.}
#'         \item{\code{demographic_stochasticity}}{Boolean for optionally choosing
#'          demographic stochasticity for the transformation.}
#'         \item{\code{density_stages}}{Array of booleans or numeric (0,1) for
#'         each stage to indicate which stages are affected by density.}
#'         \item{\code{r}}{Simulation replicate.}
#'         \item{\code{tm}}{Simulation time step.}
#'         \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'         each population at time step.}
#'         \item{\code{stage_abundance}}{Matrix of (current) abundance for each
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
#'  \item{\code{mortality}}{An optional user-defined function (optionally nested
#'  in a list with additional attributes) for applying mortality (to
#'  abundances): \code{function(params)} as per translocation.}
#'  \item{\code{dispersal}}{A list that is either length 1 or the same length as
#'  \code{stages}. If it is length 1, the same dispersal will be applied across
#'  all stages. Within each element of the list, there should be either a matrix
#'   of dispersal rates between populations (source columns to target rows) or a
#'   list of data frames of non-zero dispersal rates and indices for constructing
#'   a compact dispersal matrix, and optional changing rates over time (as per
#'   class \code{\link{DispersalGenerator}} \emph{dispersal_data} attribute).}
#'  \item{\code{dispersal_source_n_k}}{Dispersal
#'  proportion (p) density dependence via source population abundance divided by
#'  carrying capacity (n/k), where p is reduced via a linear slope (defined by
#'  two list items) from n/k <= \emph{cutoff} (p = 0) to n/k >= \emph{threshold}
#'  (aliases: \emph{dispersal_n_k_cutoff} & \emph{dispersal_n_k_threshold}).}
#'  \item{\code{dispersal_target_k}}{Dispersal rate (r) density dependence via
#'  target population carrying capacity (k), where r is reduced via a linear
#'  slope (through the origin) when k <= \emph{threshold} (alias:
#'  \emph{dispersal_k_threshold}).}
#'  \item{\code{dispersal_target_n}}{Dispersal
#'  rate (r) density dependence via target population abundance (n), where r is
#'  reduced via a linear slope (defined by two list items) from n >=
#'  \emph{threshold} to n <= \emph{cutoff} (r = 0) or vice versa (aliases:
#'  \emph{dispersal_n_threshold} & \emph{dispersal_n_cutoff}).}
#'  \item{\code{dispersal_target_n_k}}{Dispersal rate (r) density dependence via
#'  target population abundance divided by carrying capacity (n/k), where r is
#'  reduced via a linear slope (defined by two list items) from n/k >=
#'  \emph{threshold} to n/k <= \emph{cutoff} (r = 0) or vice versa.}
#'  \item{\code{season_functions}}{A list of population transformation functions
#'  (functions that change abundance across stages and compartments) the same
#'  length as \code{seasons.} The function must be in the form
#'  \code{function(params)}, where \code{params} is a list passed to the function
#'   containing:
#'       \describe{
#'          \item{\code{replicates}}{Number of replicate simulation runs
#'          (default is 1.)}
#'          \item{\code{time_steps}}{Number of simulation years. Required input.}
#'          \item{\code{seasons}}{Number of seasons per year (default is 2.)}
#'          \item{\code{populations}}{Number of populations. Required input.}
#'          \item{\code{stages}}{Number of life cycle stages. Default: 1.}
#'          \item{\code{compartments}}{Number of disease compartments (e.g., 3
#'          for a SIR model). Default: 1.}
#'          \item{\code{breeding_season_length}}{Array (matrix) of breeding season
#'           length values in days at each population cell (\code{populations}
#'           rows by \code{time_steps} columns when across time).}
#'          \item{\code{season_lengths}}{Vector of season lengths in days. Length
#'          must equal \code{seasons}.}
#'          \item{\code{mortality}}{A vector of mortality rates, one for each
#'          combination of stages and compartments. Assumed by default to be daily
#'          mortality rates. Required input.}
#'          \item{\code{mortality_unit}}{A vector indicating whether mortality
#'          rates are daily or seasonal. 1 indicates seasonal, 0 indicates daily.
#'          Default: all 0.}
#'          \item{\code{fecundity}}{A vector of fecundity rates, one for each
#'          combination of stages and compartments for which fecundity applies (see
#'          \code{fecundity_mask} below). Required input.}
#'          \item{\code{fecundity_unit}}{A vector indicating whether mortality rates
#'          are daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default:
#'          all 0.}
#'          \item{\code{fecundity_mask}}{A vector indicating which stages and
#'          compartments reproduce. Must be the same length as
#'          \code{stages * compartments}.}
#'          \item{\code{abundance_threshold}}{A quasi-extinction threshold below
#'          which a population becomes extinct.}
#'          \item{\code{demographic_stochasticity}}{Boolean for choosing demographic
#'          stochasticity for transition, dispersal, and/or other processes
#'          (default is TRUE).}
#'          \item(\code{transmission}}{A vector of transmission rates, one for each
#'          combination of stages and compartments. Assumed by default to be daily
#'          transmission rates. Required input.}
#'          \item{\code{transmission_unit}}{A vector indicating whether mortality
#'          rates are daily or seasonal. 1 indicates seasonal, 0 indicates daily.
#'          Default: all 0.}
#'          \item{\code{recovery}}{A vector of recovery rates, one for each
#'          combination of stages and compartment for which recovery applies (see
#'          \code{recovery_mask} below.)}
#'          \item{\code{recovery_unit}}{A vector indicating whether mortality rates
#'           are daily or seasonal. 1 indicates seasonal, 0 indicates daily. Default:
#'           all 0.}
#'          \item{\code{recovery_mask}}{A vector indicating which compartments are
#'          subject to recovery (i.e., infected classes that can recover.) Must be
#'          the same length as \code{compartments}.}
#'          \item{\code{r}}{Simulation replicate.} \item{\code{tm}}{Simulation time
#'          step.}
#'          \item{\code{carrying_capacity}}{Array of carrying capacity values for
#'          each population at time step.}
#'          \item{\code{stage_abundance}}{Matrix of abundance for each combination of
#'           stage and compartment (rows) and population (columns) at time step.}
#'          \item{\code{occupied_indices}}{Array of indices for populations occupied
#'          at time step.}
#'          \item{\code{simulator}}{\code{\link{SimulatorReference}} object
#'          with dynamically accessible \emph{attached} and \emph{results} lists.}
#'          \item{\code{additional attributes}}{Additional attributes when the
#'          transformation is optionally nested in a list.}}
#'       and returns a transformed stage abundance matrix.
#'  }
#'  \item{\code{simulation_order}}{A list the same length as \code{seasons}.
#'  Each element in the list is a vector of named simulation processes in the
#'  desired order. Processes must be one of "transition", "translocation",
#'  "harvest", "mortality", "dispersal", "season_functions", or "results."
#'  "season_functions" will be matched to the appropriate season (i.e., if
#'  "season_functions" appears in element 1 of the list, season_functions[[1]]
#'  will be called.) Required input.
#'  \item{\code{additional transformation functions}}{Additional user-defined
#'  abundance transformation functions (optionally nested in lists with
#'  additional attributes) are utilised when listed in \emph{simulation_order}
#'  (function as per translocation).}
#'  \item{\code{results_selection}}{List of results selection from: "abundance"
#'  (default), "ema", "extirpation", "extinction_location", "harvested",
#'  "occupancy"; "summarize" (default) or "replicate".}
#'  \item{\code{result_stages}}{Array of booleans or numeric (0, 1, 2, ...) for
#'   each stage to indicate which stages are included/combined (each unique digit
#'    > 0; optionally named) in the results (default is 1 for all stages).}
#'  \item{\code{result_compartments}}{Array of booleans or numeric for each
#'  compartment to indicate which stages are included/combined (each unique digit
#'   > 0; optionally named) in the results (default is 1 for all compartments).}
#'}
#'@return Selected simulation results as a nested list summarized (mean, sd,
#'  min, max) across multiple replicates (default), or 2-3D arrays including
#'  results for each replicate:
#' \describe{
#'   \item{\code{abundance}}{Matrix or 3D array of simulation abundance:
#'   \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'   deep).}
#'     \item{\code{abundance_stages}}{List of matrices or 3D arrays of simulation
#'     abundance for unique stage-compartment combinations when present: each
#'     \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'     deep).}
#'     \item{\code{all$abundance}}{Array or matrix of total abundance across
#'     populations: \emph{time_steps} (rows by \emph{replicates} columns).}
#'     \item{\code{all$abundance_stages}}{List of arrays or matrices of total
#'     abundance across populations for unique stage-compartment combinations when
#'      present: each \emph{time_steps} (rows by \emph{replicates} columns).}
#'     \item{\code{all$ema}}{Array of expected minimum abundance at each time step
#'      (averaged across replicates).}
#'     \item{\code{extirpation}}{Array or matrix of extirpation times:
#'     \emph{populations} (rows by \emph{replicates} columns).}
#'     \item{\code{all$extirpation}}{Array of extirpation time across populations
#'      for each replicate.}
#'     \item{\code{all$extinction_location}}{The weighted centroid of cells
#'     occupied in the time-step prior to the extirpation of all populations
#'     (if it occurred) for each replicate.}
#'     \item{\code{harvested}}{Matrix or 3D array of individuals harvested:
#'     \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'     deep).}
#'     \item{\code{harvested_stages}}{List of matrices or 3D arrays of individuals
#'      harvested for unique stage-compartment combinations when present: each
#'      \emph{populations} rows by \emph{time_steps} columns (by \emph{replicates}
#'       deep).}
#'     \item{\code{all$harvested}}{Array or matrix of individuals harvested across
#'      populations: \emph{time_steps} (rows by \emph{replicates} columns).}
#'     \item{\code{all$harvested_stages}}{List of arrays or matrices of individuals
#'      harvested across populations for unique stage-compartment combinations when
#'       present: each \emph{time_steps} (rows by \emph{replicates} columns).}
#'     \item{\code{all$occupancy}}{Array or matrix of the number of populations
#'     occupied at each time-step: \emph{time_steps} (rows by \emph{replicates}
#'     columns).}
#'     \item{\code{additional results}}{Additional results may be attached via
#'     user-defined functions (using \code{params$simulator$results}).}
#' }
#'
#' @import poems
#' @import cli
#' @export disease_simulator

disease_simulator <- function(inputs) {

  # Stop if minimal inputs are not present
  if (is.null(inputs[["time_steps"]]) || is.null(inputs[["populations"]]) ||
      is.null(inputs[["initial_abundance"]]) ||
      is.null(inputs[["carrying_capacity"]]) ||
      is.null(inputs[["mortality"]]) || is.null(inputs[["transmission"]]) ||
      is.null(inputs[["fecundity"]]) || is.null(inputs[["simulation_order"]])) {
    incomplete_inputs <- if (is.null(inputs[["time_steps"]]))
      "time_steps"
    incomplete_inputs <- c(incomplete_inputs, if (is.null(inputs[["populations"]]))
        "populations")
    incomplete_inputs <- c(incomplete_inputs, if (is.null(inputs[["initial_abundance"]]))
        "initial_abundance")
    incomplete_inputs <- c(incomplete_inputs, if (is.null(inputs[["carrying_capacity"]]))
      "carrying_capacity")
    incomplete_inputs <- c(incomplete_inputs, if (is.null(inputs[["fecundity"]]))
        "fecundity")
    cli_abort(c("Minimal inputs required to run simulation are missing.",
                "x" = "Your input does not include {incomplete_inputs}."))
  }

  ## Unpack inputs and calculate/initialize re-usable variables ##

  # General simulation settings
  if (!is.null(inputs[["random_seed"]])) {
    set.seed(inputs[["random_seed"]])
  }
  replicates <- ifelse(is.null(inputs[["replicates"]]), 1, inputs[["replicates"]])
  time_steps <- inputs[["time_steps"]]
  years_per_step <- 1 # because some poems functions require this
  seasons <- ifelse(is.null(inputs[["seasons"]]), 2, inputs[["seasons"]])
  populations <- inputs[["populations"]]
  stages <- ifelse(is.null(inputs[["stages"]]), 1, inputs[["stages"]])
  compartments <- ifelse(is.null(inputs[["compartments"]]), 1, inputs[["compartments"]])

  # Initial abundance for each stage and population
  if (inherits(inputs, "DiseaseModel") &&
      !is.null(inputs[["region"]]) &&
      inputs[["region"]][["use_raster"]] &&
      inherits(class(inputs[["initial_abundance"]]),
               c("RasterLayer", "RasterStack", "RasterBrick"))) {
    initial_abundance <- as.matrix(
      inputs[["initial_abundance"]][inputs[["region"]][["region_indices"]]]
    )
  } else if (is.vector(inputs[["initial_abundance"]]) &&
             is.numeric(inputs[["initial_abundance"]])) {
    initial_abundance <- matrix(inputs[["initial_abundance"]], nrow = 1)
  } else if (is.array(inputs[["initial_abundance"]])) {
    initial_abundance <- inputs[["initial_abundance"]]
  } else {
    initial_abundance <- inputs[["initial_abundance"]]
    cli_abort(
      c(
        "{.var initial_abundance} must be a numeric vector, raster, matrix, or
                array.",
        "x" = "{.var initial_abundance} is a {.obj_type_friendly {initial_abundance}}."
      )
    )
  }

  segments <- stages*compartments

  if (nrow(initial_abundance) != segments) {
    cli_abort(c("{.var initial_abundance} has {nrow(initial_abundance)} row{?s}.",
                "x" = "There should be {segments} row{?s} (compartments x stages)."))
  }
  if (ncol(initial_abundance) != populations) {
    cli_abort(c("{.var initial_abundance} has {ncol(initial_abundance)} column{?s}.",
                "x" = "There should be {populations} column{?s} ({.var populations})."))
  }
  population_abundance <- .colSums(initial_abundance, m = segments, n = populations)

  # Check carrying capacity

}