#' Nested functions for initializing, calculating and collecting disease
#' simulator results.
#'
#' Modular functions for the disease simulator for initializing, calculating and
#' collecting simulator results.
#'
#' @param replicates Number of replicate simulation runs.
#' @param time_steps Number of simulation time steps.
#' @param seasons Number of seasons per time step.
#' @param coordinates Data frame (or matrix) of X-Y population coordinates.
#' @param initial_abundance Matrix of initial abundances at each combination of
#' stage and compartment (in rows) for each population (in columns).
#' @param results_selection List of results selection from: "abundance"
#' (default), "ema", "extirpation", "extinction_location", "harvested",
#' "occupancy"; "summarize" (default) or "replicate".
#' @param result_stages Array of booleans or numeric (0, 1, 2, ...) for each
#' stage to indicate which stages are included/combined (each unique digit > 0;
#' optionally named) in the results (default is 1 for all stages).
#' @param result_compartments Array of booleans or numeric for each
#'  compartment to indicate which stages are included/combined (each unique
#'  digit \> 0; optionally named) in the results (default is 1 for all
#'  compartments).
#' @return List of result functions:
#'   \describe{
#'     \item{\code{initialize_attributes = function())}}{Constructs and returns
#'     an initialized nested list for the selected result attributes.}
#'     \item{\code{initialize_replicate = function(results)}}{Initializes and
#'     returns nested result attributes at the start of each replicate.}
#'     \item{\code{calculate_at_season = function(r, tm, season,
#'     segment_abundance, harvested, results)}}{Appends and calculates
#'     (non-NULL) results and returns nested result attributes at the end of
#'     each season within time step (tm) within replicate (r).}
#'     \item{\code{calculate_at_timestep = function(r, tm, segment_abundance,
#'     harvested, results)}}{Appends and calculates (non-NULL) results and
#'     returns nested result attributes at the end of each time step (tm) within
#'      replicate (r).}
#'     \item{\code{finalize_attributes = function(results)}}{Finalizes result
#'     calculations at the end of the simulation.}
#'   }
#' @export disease_results

disease_results <- function(replicates, time_steps, seasons, coordinates,
                            initial_abundance, results_selection = NULL,
                            result_stages = NULL, result_compartments = NULL) {

  # Set defaults when NULL
  if (is.null(results_selection)) {
    results_selection <- c("abundance", "summarize")
  }


}
