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

disease_results <- function(replicates,
                            time_steps,
                            seasons,
                            coordinates,
                            initial_abundance,
                            results_selection = NULL,
                            result_stages = NULL,
                            result_compartments = NULL) {

  # Set defaults when NULL
  if (is.null(results_selection)) {
    results_selection <- c("abundance", "summarize")
  }
  separated_stage_indices <- NULL
  if (!is.null(result_stages)) {
    result_stage_indices <- which(result_stages > 0)
    if (length(unique(result_stages[result_stage_indices])) > 1) {
      # separated stage combinations for abundance and harvest results
      separated_stage_indices <-
        apply(as.array(unique(result_stages[result_stage_indices])), 1, function(s)
          which(result_stages == s))
      names(separated_stage_indices) <- lapply(separated_stage_indices, function(s)
        paste(unique(names(s)), collapse = "."))
    }
  }
  result_stage_length <- length(result_stage_indices)
  separated_compartment_indices <- NULL
  if (!is.null(result_compartments)) {
    result_compartment_indices <- which(result_compartments > 0)
    if (length(unique(result_compartments[result_compartment_indices])) > 1) {
      # separated compartment combinations for abundance and harvest results
      separated_compartment_indices <-
        apply(as.array(unique(result_compartments[result_compartment_indices])), 1, function(s)
          which(result_compartments == s))
      names(separated_compartment_indices) <- lapply(separated_compartment_indices, function(s)
        paste(unique(names(s)), collapse = "."))
    }
  }
  result_compartment_length <- length(result_compartment_indices)
  if (is.null(result_stages) | is.null(result_compartments)) {
    segment_indices <- 1:nrow(initial_abundance)
  } else {
    selected_indices <- expand.grid(result_stage_indices,
                                    result_compartment_indices)
    segment_indices <- sort((selected_indices$Var1 - 1) * result_compartment_length +
                              selected_indices$Var2)
  }

  # Initialize abundance count
  populations <- ncol(initial_abundance)
  initial_abundance_count <- as.array(
    .colSums(initial_abundance[segment_indices,],
             m = result_stage_length*result_compartment_length,
             n = populations)
  )

  # Maintain initial abundance count for EMA calculations
  if ("ema" %in% results_selection) {
    initial_abundance_count_min <- sum(initial_abundance_count)
  }

  initialize_attributes <- function() {

    # Initialize lists for results for each population and for all populations
    results <- list()
    results$all <- list()

    # Abundance results
    if ("abundance" %in% results_selection) {
      results$abundance <- array(0, c(populations, time_steps, seasons))
      results$all$abundance <- array(0, c(time_steps, seasons))
      if (replicates > 1) {
        if ("replicate" %in% results_selection) {
          results$abundance <- array(0, c(populations, time_steps, seasons,
                                          replicates))
          results$all$abundance <- array(0, c(time_steps, seasons, replicates))
        } else { # summarize (default)
          results$abundance <- list(mean = results$abundance,
                                    sd = results$abundance,
                                    min = results$abundance,
                                    max = results$abundance)
          results$all$abundance <- list(mean = results$all$abundance,
                                        sd = results$all$abundance,
                                        min = results$all$abundance,
                                        max = results$all$abundance)
        }
      }
      if (length(separated_stage_indices) && 
      length(separated_compartment_indices)) {
        results$abundance_segments <- lapply(segment_indices,
                                             function(s) results$abundance)
        results$all$abundance_segments <- lapply(segment_indices,
                                                 function(s) results$all$abundance)
      } else if (length(separated_stage_indices)) {
        results$abundance_stages <- lapply(separated_stage_indices,
                                           function(s) results$abundance)
        results$all$abundance_stages <- lapply(separated_stage_indices,
                                               function(s) results$all$abundance)
      } else if (length(separated_compartment_indices)) {
        results$abundance_compartments <- lapply(separated_compartment_indices,
                                                 function(s) results$abundance)
        results$all$abundance_compartments <- lapply(separated_compartment_indices,
                                                     function(s) results$all$abundance)
      }
    }

    # Expected minimum abundance (EMA) results
    if ("ema" %in% results_selection) {
      results$all$ema <- array(0, time_steps, seasons)
    }

    # Extirpation results
    if ("extirpation" %in% results_selection) {
      results$extirpation <- array(as.numeric(NA), populations)
      results$extirpation[which(initial_abundance_count == 0)] <- 0
      results$all$extirpation <- array(as.numeric(NA), replicates)
      results$all$extirpation[all(initial_abundance_count == 0)] <- 0
      if (replicates > 1) {
        results$extirpation <- array(results$extirpation, c(populations, replicates))
      }
    }

    # Extinction location results
    if ("extinction_location" %in% results_selection) {
      results$all$extinction_location <- array(as.numeric(NA), c(replicates, 2))
      colnames(results$all$extinction_location) <- c("x", "y")
      results$last_occupied_abundance_count <- initial_abundance_count
    }

    # Harvested results
    if ("harvested" %in% results_selection) {
      results$harvested <- array(0, c(populations, time_steps, seasons))
      results$all$harvested <- array(0, c(time_steps, seasons))
      if (replicates > 1) {
        if ("replicate" %in% results_selection) {
          results$harvested <- array(0, c(populations, time_steps, seasons,
                                          replicates))
          results$all$harvested <- array(0, c(time_steps, seasons, replicates))
        } else { # summarize (default)
          results$harvested <- list(mean = results$harvested,
                                    sd = results$harvested,
                                    min = results$harvested,
                                    max = results$harvested)
          results$all$harvested <- list(mean = results$all$harvested,
                                        sd = results$all$harvested,
                                        min = results$all$harvested,
                                        max = results$all$harvested)
        }
      }
      if (length(separated_stage_indices) && length(separated_compartment_indices)) {
        results$harvested_segments <- lapply(segment_indices,
                                             function(s) results$harvested)
        results$all$harvested_segments <- lapply(segment_indices,
                                                 function(s) results$all$harvested)
      } else if (length(separated_stage_indices)) {
        results$harvested_stages <- lapply(separated_stage_indices,
                                           function(s) results$harvested)
        results$all$harvested_stages <- lapply(separated_stage_indices,
                                               function(s) results$all$harvested)
      } else if (length(separated_compartment_indices)) {
        results$harvested_compartments <- lapply(separated_compartment_indices,
                                                 function(s) results$harvested)
        results$all$harvested_compartments <- lapply(separated_compartment_indices,
                                                     function(s) results$all$harvested)
      }
    }

    # Occupancy results
    if ("occupancy" %in% results_selection) {
      results$all$occupancy <- array(0, c(time_steps, seasons))
      if (replicates > 1) {
        if ("replicate" %in% results_selection) {
          results$all$occupancy <- array(0, c(time_steps, seasons, replicates))
        } else { # summarize (default)
          results$all$occupancy <- list(mean = results$all$occupancy,
                                        sd = results$all$occupancy,
                                        min = results$all$occupancy,
                                        max = results$all$occupancy)
        }
      }
    }

    return(results)
  } # End initialize attributes

  # Create a nested function for (re-)initializing the results attributes at
  # the beginning of each replicate
  initialize_replicate <- function(results) {
    if ("ema" %in% results_selection) {
      results$abundance_count_min <- initial_abundance_count_min
    }
    if ("extinction_location" %in% results_selection) {
      results$last_occupied_abundance_count <- initial_abundance_count
    }
    return(results)
  }

  return(list(initialize_attributes = initialize_attributes,
              initialize_replicate = initialize_replicate,
              calculate_at_timestep = calculate_at_timestep,
              calculate_at_replicate = calculate_at_replicate,
              finalize_attributes = finalize_attributes))
} # End disease results function
