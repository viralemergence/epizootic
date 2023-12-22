#' Nested functions for initializing, calculating and collecting disease
#' simulator results.
#'
#' Modular functions for the disease simulator for initializing, calculating and
#' collecting simulator results.
#'
#' @param replicates Number of replicate simulation runs.
#' @param time_steps Number of simulation time steps.
#' @param seasons Number of seasons per time step.
#' @param stages Number of life cycle stages.
#' @param compartments Number of disease compartments.
#' @param coordinates Data frame (or matrix) of X-Y population coordinates.
#' @param initial_abundance Matrix of initial abundances at each combination of
#' stage and compartment (in rows) for each population (in columns).
#' @param results_selection List of results selection from: "abundance"
#' (default), "ema", "extirpation", "extinction_location", "harvested",
#' "occupancy"; "summarize" (default) or "replicate". "summarize" calculates
#' mean, sd, min and max across replicates. "replicate" returns results
#' separately for each replicate.
#' @param results_breakdown A string with one of these values: "segments" (default),
#' "compartments", "stages" or "pooled." "segments" returns results for each
#' segment (stage x compartment combination.) "compartments" returns results for
#' each disease compartment. "stages" returns results for each life cycle stage.
#' "pooled" returns results that are not broken down by stage or compartment.
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
                            stages,
                            compartments,
                            coordinates,
                            initial_abundance,
                            results_selection = NULL,
                            results_breakdown = NULL) {
  # Set defaults when NULL
  if (is.null(results_selection)) {
    results_selection <- c("abundance", "summarize")
  }
  if (is.null(results_breakdown)) {
    results_breakdown <- c("segments")
  }

  # Initialize abundance count
  populations <- ncol(initial_abundance)
  initial_abundance_count <- as.array(.colSums(initial_abundance,
                                               m = stages * compartments,
                                               n = populations))

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
          results$abundance <-
            array(0, c(populations, time_steps, seasons, replicates))
          results$all$abundance <-
            array(0, c(time_steps, seasons, replicates))
        } else {
          # summarize (default)
          results$abundance <- list(
            mean = results$abundance,
            sd = results$abundance,
            min = results$abundance,
            max = results$abundance
          )
          results$all$abundance <- list(
            mean = results$all$abundance,
            sd = results$all$abundance,
            min = results$all$abundance,
            max = results$all$abundance
          )
        }
      }
      if (results_breakdown == "segments") {
        results$abundance_segments <- lapply(1:(stages * compartments),
                                             function(s)
                                               results$abundance)
        results$all$abundance_segments <- lapply(1:(stages * compartments),
                                                 function(s)
                                                   results$all$abundance)
        names(results$abundance_segments) <- expand.grid(1:stages, 1:compartments) |>
          pmap_chr(function(Var1, Var2) {
            paste("stage", Var1, "compartment", Var2, sep = "_")
          })
        names(results$all$abundance_segments) <- names(results$abundance_segments)
      } else if (results_breakdown == "stages") {
        results$abundance_stages <- lapply(1:stages,
                                           function(s)
                                             results$abundance)
        results$all$abundance_stages <- lapply(1:stages,
                                               function(s)
                                                 results$all$abundance)
        names(results$abundance_stages) <- map_chr(1:stages, function(s) {
          paste("stage", s, sep = "_")
        })
        names(results$all$abundance_stages) <- names(results$abundance_stages)
      } else if (results_breakdown == "compartments") {
        results$abundance_compartments <- lapply(1:compartments,
                                                 function(s)
                                                   results$abundance)
        results$all$abundance_compartments <- lapply(1:compartments,
                                                     function(s)
                                                       results$all$abundance)
        names(results$abundance_compartments) <- map_chr(1:compartments, function(s) {
          paste("compartment", s, sep = "_")
        })
        names(results$all$abundance_compartments) <- names(results$abundance_compartments)
      }
    }

    # Expected minimum abundance (EMA) results
    if ("ema" %in% results_selection) {
      results$all$ema <- array(0, c(time_steps, seasons))
    }

    # Extirpation results
    if ("extirpation" %in% results_selection) {
      results$extirpation <- array(as.numeric(NA), populations)
      results$extirpation[which(initial_abundance_count == 0)] <- 0
      results$all$extirpation <- array(as.numeric(NA), replicates)
      results$all$extirpation[all(initial_abundance_count == 0)] <-
        0
      if (replicates > 1) {
        results$extirpation <- array(results$extirpation,
                                     c(populations, replicates))
      }
    }

    # Extinction location results
    if ("extinction_location" %in% results_selection) {
      results$all$extinction_location <-
        array(as.numeric(NA), c(replicates, 2))
      colnames(results$all$extinction_location) <- c("x", "y")
      results$last_occupied_abundance_count <-
        initial_abundance_count
    }

    # Harvested results
    if ("harvested" %in% results_selection) {
      results$harvested <- array(0, c(populations, time_steps, seasons))
      results$all$harvested <- array(0, c(time_steps, seasons))
      if (replicates > 1) {
        if ("replicate" %in% results_selection) {
          results$harvested <- array(0, c(populations, time_steps, seasons,
                                          replicates))
          results$all$harvested <-
            array(0, c(time_steps, seasons, replicates))
        } else {
          # summarize (default)
          results$harvested <- list(
            mean = results$harvested,
            sd = results$harvested,
            min = results$harvested,
            max = results$harvested
          )
          results$all$harvested <- list(
            mean = results$all$harvested,
            sd = results$all$harvested,
            min = results$all$harvested,
            max = results$all$harvested
          )
        }
      }
      if (results_breakdown == "segments") {
        results$harvested_segments <- lapply(1:(stages * compartments),
                                             function(s)
                                               results$harvested)
        names(results$harvested_segments) <-
          expand.grid(1:stages, 1:compartments) |>
          pmap_chr(function(Var1, Var2) {
            paste("stage", Var1, "compartment", Var2, sep = "_")
          })
        results$all$harvested_segments <- lapply(1:(stages * compartments),
                                                 function(s)
                                                   results$all$harvested)
        names(results$all$harvested_segments) <-
          names(results$harvested_segments)
      } else if (results_breakdown == "stages") {
        results$harvested_stages <- lapply(1:stages,
                                           function(s)
                                             results$harvested)
        names(results$harvested_stages) <- map_chr(1:stages, function(s) {
          paste("stage", s, sep = "_")
        })
        results$all$harvested_stages <- lapply(1:stages,
                                               function(s)
                                                 results$all$harvested)
        names(results$all$harvested_stages) <-
          names(results$harvested_stages)
      } else if (results_breakdown == "compartments") {
        results$harvested_compartments <- lapply(1:compartments,
                                                 function(s)
                                                   results$harvested)
        names(results$harvested_compartments) <- map_chr(1:compartments, function(s) {
          paste("compartment", s, sep = "_")
        })
        results$all$harvested_compartments <- lapply(1:compartments,
                                                     function(s)
                                                       results$all$harvested)
        names(results$all$harvested_compartments) <-
          names(results$harvested_compartments)
      }
    }

    # Occupancy results
    if ("occupancy" %in% results_selection) {
      results$all$occupancy <- array(0, c(time_steps, seasons))
      if (replicates > 1) {
        if ("replicate" %in% results_selection) {
          results$all$occupancy <-
            array(0, c(time_steps, seasons, replicates))
        } else {
          # summarize (default)
          results$all$occupancy <- list(
            mean = results$all$occupancy,
            sd = results$all$occupancy,
            min = results$all$occupancy,
            max = results$all$occupancy
          )
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

  ## Create a nested function for calculating/collecting (non-NULL) results at season
  calculate_at_season <-
    function(r,
             tm,
             season,
             segment_abundance,
             harvested,
             results) {
      if (!is.null(segment_abundance)) {
        # results calculated from abundance

        # Select abundance counts for included stages
        abundance_count <- as.array(.colSums(
          segment_abundance,
          m = stages * compartments,
          n = populations
        ))
        all_abundance_count <- sum(abundance_count)

        if (results_breakdown == "segments") {
          abundance_segments_count <-
            lapply(1:(stages * compartments), function(s) {
              as.array(.colSums(segment_abundance[s,], m = length(s), n = populations))
            })
          all_abundance_segments_count <-
            lapply(abundance_segments_count, sum)
        } else if (results_breakdown == "stages") {
          abundance_stages_count <- lapply(1:stages, function(s) {
            stage_indices <- seq(s, nrow(segment_abundance), by = stages)
            as.array(.colSums(
              segment_abundance[stage_indices,],
              m = length(stage_indices),
              n = populations
            ))
          })
          all_abundance_stages_count <-
            lapply(abundance_stages_count, sum)
        } else if (results_breakdown == "compartments") {
          abundance_compartments_count <- lapply(1:compartments, function(s) {
            compartment_indices <- seq(stages * s - (stages - 1), stages * s, 1)
            as.array(.colSums(
              segment_abundance[compartment_indices,],
              m = length(compartment_indices),
              n = populations
            ))
          })
          all_abundance_compartments_count <-
            lapply(abundance_compartments_count, sum)
        }

        # Abundance results
        if ("abundance" %in% results_selection) {
          if (replicates > 1) {
            if ("replicate" %in% results_selection) {
              results$abundance[, tm, season, r] <- abundance_count
              results$all$abundance[tm, season, r] <-
                all_abundance_count
            } else {
              # summarize (default)
              previous_abundance_mean <-
                results$abundance$mean[, tm, season]
              results$abundance$mean[, tm, season] <-
                previous_abundance_mean +
                (abundance_count - previous_abundance_mean) / r
              previous_abundance_var_by_r <-
                results$abundance$sd[, tm, season]
              results$abundance$sd[, tm, season] <-
                previous_abundance_var_by_r +
                (abundance_count - previous_abundance_mean) *
                (abundance_count - results$abundance$mean[, tm, season])

              if (r > 1) {
                results$abundance$min[, tm, season] <- pmin(results$abundance$min[, tm, season],
                                                            abundance_count)
                results$abundance$max[, tm, season] <- pmax(results$abundance$max[, tm, season],
                                                            abundance_count)
              } else {
                results$abundance$min[, tm, season] <- abundance_count
                results$abundance$max[, tm, season] <- abundance_count
              }
              previous_all_abundance_mean <-
                results$all$abundance$mean[tm, season]
              results$all$abundance$mean[tm, season] <-
                previous_all_abundance_mean +
                (all_abundance_count - previous_all_abundance_mean) / r
              previous_all_abundance_var_by_r <-
                results$all$abundance$sd[tm, season]
              results$all$abundance$sd[tm, season] <-
                previous_all_abundance_var_by_r +
                (all_abundance_count - previous_all_abundance_mean) *
                (all_abundance_count - results$all$abundance$mean[tm, season])
              results$all$abundance$min[tm, season] <- ifelse(
                r > 1,
                min(results$all$abundance$min[tm, season], all_abundance_count),
                all_abundance_count
              )
              results$all$abundance$max[tm, season] <- ifelse(
                r > 1,
                max(results$all$abundance$max[tm, season], all_abundance_count),
                all_abundance_count
              )
            }
          } else {
            # single replicate
            results$abundance[, tm, season] <- abundance_count
            results$all$abundance[tm, season] <- all_abundance_count
          }
          # Update abundance and all abundance for each segment
          if (results_breakdown == "segments") {
            segment_names <- names(results$abundance_segments)
            for (i in 1:(stages * compartments)) {
              if (replicates > 1) {
                if ("replicate" %in% results_selection) {
                  results$abundance_segments[[i]][, tm, season, r] <-
                    abundance_segments_count[[i]]
                  results$all$abundance_segments[[i]][tm, season, r] <-
                    all_abundance_segments_count[[i]]
                } else {
                  # Update mean and standard deviation for each segment
                  previous_abundance_segments_mean <-
                    results$abundance_segments[[i]]$mean[, tm, season]
                  results$abundance_segments[[i]]$mean[, tm, season] <-
                    previous_abundance_segments_mean +
                    (abundance_segments_count[[i]] - previous_abundance_segments_mean) / r
                  previous_abundance_segments_var_by_r <-
                    results$abundance_segments[[i]]$sd[, tm, season]
                  results$abundance_segments[[i]]$sd[, tm, season] <-
                    previous_abundance_segments_var_by_r +
                    (abundance_segments_count[[i]] - previous_abundance_segments_mean) *
                    (abundance_segments_count[[i]] - results$abundance_segments[[i]]$mean[, tm, season])
                  # Update min and max for each segment
                  if (r == 1) {
                    results$abundance_segments[[i]]$min[, tm, season] <-
                      abundance_segments_count[[i]]
                    results$abundance_segments[[i]]$max[, tm, season] <-
                      abundance_segments_count[[i]]
                  } else {
                    results$abundance_segments[[i]]$min[, tm, season] <- pmin(results$abundance_segments[[i]]$min[, tm, season],
                                                                              abundance_segments_count[[i]])
                    results$abundance_segments[[i]]$max[, tm, season] <-
                      pmax(results$abundance_segments[[i]]$max[, tm, season],
                           abundance_segments_count[[i]])
                  }

                  previous_all_abundance_segments_mean <-
                    results$all$abundance_segments[[i]]$mean[tm, season]
                  results$all$abundance_segments[[i]]$mean[tm, season] <-
                    previous_all_abundance_segments_mean +
                    (
                      all_abundance_segments_count[[i]] - previous_all_abundance_segments_mean
                    ) / r
                  previous_all_abundance_segments_var_by_r <-
                    results$all$abundance_segments[[i]]$sd[tm, season]
                  results$all$abundance_segments[[i]]$sd[tm, season] <-
                    previous_all_abundance_segments_var_by_r +
                    (
                      all_abundance_segments_count[[i]] - previous_all_abundance_segments_mean
                    ) *
                    (
                      all_abundance_segments_count[[i]] - results$all$abundance_segments[[i]]$mean[tm, season]
                    )

                  results$all$abundance_segments[[i]]$min[tm, season] <-
                    ifelse(
                      r > 1,
                      min(
                        results$all$abundance_segments[[i]]$min[tm, season],
                        all_abundance_segments_count[[i]]
                      ),
                      all_abundance_segments_count[[i]]
                    )
                  results$all$abundance_segments[[i]]$max[tm, season] <-
                    ifelse(
                      r > 1,
                      max(
                        results$all$abundance_segments[[i]]$max[tm, season],
                        all_abundance_segments_count[[i]]
                      ),
                      all_abundance_segments_count[[i]]
                    )
                }
              } else {
                # single replicate
                results$abundance_segments[[i]][, tm, season] <-
                  abundance_segments_count[[i]]
                results$all$abundance_segments[[i]][tm, season] <-
                  all_abundance_segments_count[[i]]
              }
            }
            names(results$all$abundance_segments) <- segment_names
            names(results$abundance_segments) <- segment_names
          }

          # Update abundance and all abundance for each stage
          if (results_breakdown == "stages") {
            stage_names <- names(results$abundance_stages)
            for (i in 1:stages) {
              if (replicates > 1) {
                if ("replicate" %in% results_selection) {
                  results$abundance_stages[[i]][, tm, season, r] <- abundance_stages_count[[i]]
                  results$all$abundance_stages[[i]][tm, season, r] <- all_abundance_stages_count[[i]]
                } else {
                  # Update mean and standard deviation for each stage
                  previous_abundance_stages_mean <- results$abundance_stages[[i]]$mean[, tm, season]
                  results$abundance_stages[[i]]$mean[, tm, season] <- previous_abundance_stages_mean +
                    (abundance_stages_count[[i]] - previous_abundance_stages_mean) / r
                  previous_abundance_stages_var_by_r <- results$abundance_stages[[i]]$sd[, tm, season]
                  results$abundance_stages[[i]]$sd[, tm, season] <- previous_abundance_stages_var_by_r +
                    (abundance_stages_count[[i]] - previous_abundance_stages_mean) *
                    (abundance_stages_count[[i]] - results$abundance_stages[[i]]$mean[, tm, season])
                  # Update min and max for each stage
                  if (r == 1) {
                    results$abundance_stages[[i]]$min[, tm, season] <- abundance_stages_count[[i]]
                    results$abundance_stages[[i]]$max[, tm, season] <- abundance_stages_count[[i]]
                  } else {
                    results$abundance_stages[[i]]$min[, tm, season] <- pmin(results$abundance_stages[[i]]$min[, tm, season],
                                                                            abundance_stages_count[[i]])
                    results$abundance_stages[[i]]$max[, tm, season] <- pmax(results$abundance_stages[[i]]$max[, tm, season],
                                                                            abundance_stages_count[[i]])
                  }

                  previous_all_abundance_stages_mean <- results$all$abundance_stages[[i]]$mean[tm, season]
                  results$all$abundance_stages[[i]]$mean[tm, season] <- previous_all_abundance_stages_mean +
                    (all_abundance_stages_count[[i]] - previous_all_abundance_stages_mean) / r
                  previous_all_abundance_stages_var_by_r <- results$all$abundance_stages[[i]]$sd[tm, season]
                  results$all$abundance_stages[[i]]$sd[tm, season] <- previous_all_abundance_stages_var_by_r +
                    (all_abundance_stages_count[[i]] - previous_all_abundance_stages_mean) *
                    (all_abundance_stages_count[[i]] - results$all$abundance_stages[[i]]$mean[tm, season])

                  results$all$abundance_stages[[i]]$min[tm, season] <-
                    ifelse(
                      r > 1,
                      min(results$all$abundance_stages[[i]]$min[tm, season], all_abundance_stages_count[[i]]),
                      all_abundance_stages_count[[i]]
                    )
                  results$all$abundance_stages[[i]]$max[tm, season] <-
                    ifelse(
                      r > 1,
                      max(results$all$abundance_stages[[i]]$max[tm, season], all_abundance_stages_count[[i]]),
                      all_abundance_stages_count[[i]]
                    )
                }
              } else {
                # single replicate
                results$abundance_stages[[i]][, tm, season] <- abundance_stages_count[[i]]
                results$all$abundance_stages[[i]][tm, season] <- all_abundance_stages_count[[i]]
              }
            }
            names(results$all$abundance_stages) <- stage_names
            names(results$abundance_stages) <- stage_names
          }

          # Update abundance and all abundance for each compartment
          if (results_breakdown == "compartments") {
            compartment_names <- names(results$abundance_compartments)
            for (i in 1:compartments) {
              if (replicates > 1) {
                if ("replicate" %in% results_selection) {
                  results$abundance_compartments[[i]][, tm, season, r] <-
                    abundance_compartments_count[[i]]
                  results$all$abundance_compartments[[i]][tm, season, r] <-
                    all_abundance_compartments_count[[i]]
                } else {
                  # Update mean and standard deviation for each compartment
                  previous_abundance_compartments_mean <-
                    results$abundance_compartments[[i]]$mean[, tm, season]
                  results$abundance_compartments[[i]]$mean[, tm, season] <-
                    previous_abundance_compartments_mean +
                    (abundance_compartments_count[[i]] - previous_abundance_compartments_mean) / r
                  previous_abundance_compartments_var_by_r <-
                    results$abundance_compartments[[i]]$sd[, tm, season]
                  results$abundance_compartments[[i]]$sd[, tm, season] <-
                    previous_abundance_compartments_var_by_r +
                    (abundance_compartments_count[[i]] - previous_abundance_compartments_mean) *
                    (abundance_compartments_count[[i]] - results$abundance_compartments[[i]]$mean[, tm, season])
                  # Update min and max for each compartment
                  if (r == 1) {
                    results$abundance_compartments[[i]]$min[, tm, season] <-
                      abundance_compartments_count[[i]]
                    results$abundance_compartments[[i]]$max[, tm, season] <-
                      abundance_compartments_count[[i]]
                  } else {
                    results$abundance_compartments[[i]]$min[, tm, season] <- pmin(
                      results$abundance_compartments[[i]]$min[, tm, season],
                      abundance_compartments_count[[i]]
                    )
                    results$abundance_compartments[[i]]$max[, tm, season] <-
                      pmax(
                        results$abundance_compartments[[i]]$max[, tm, season],
                        abundance_compartments_count[[i]]
                      )
                  }

                  previous_all_abundance_compartments_mean <-
                    results$all$abundance_compartments[[i]]$mean[tm, season]
                  results$all$abundance_compartments[[i]]$mean[tm, season] <-
                    previous_all_abundance_compartments_mean +
                    (all_abundance_compartments_count[[i]] - previous_all_abundance_compartments_mean) / r
                  previous_all_abundance_compartments_var_by_r <-
                    results$all$abundance_compartments[[i]]$sd[tm, season]
                  results$all$abundance_compartments[[i]]$sd[tm, season] <-
                    previous_all_abundance_compartments_var_by_r +
                    (all_abundance_compartments_count[[i]] - previous_all_abundance_compartments_mean) *
                    (all_abundance_compartments_count[[i]] - results$all$abundance_compartments[[i]]$mean[tm, season])

                  results$all$abundance_compartments[[i]]$min[tm, season] <-
                    ifelse(
                      r > 1,
                      min(
                        results$all$abundance_compartments[[i]]$min[tm, season],
                        all_abundance_compartments_count[[i]]
                      ),
                      all_abundance_compartments_count[[i]]
                    )
                  results$all$abundance_compartments[[i]]$max[tm, season] <-
                    ifelse(
                      r > 1,
                      max(
                        results$all$abundance_compartments[[i]]$max[tm, season],
                        all_abundance_compartments_count[[i]]
                      ),
                      all_abundance_compartments_count[[i]]
                    )
                }
              } else {
                # single replicate
                results$abundance_compartments[[i]][, tm, season] <-
                  abundance_compartments_count[[i]]
                results$all$abundance_compartments[[i]][tm, season] <-
                  all_abundance_compartments_count[[i]]
              }
            }
            names(results$all$abundance_compartments) <- compartment_names
            names(results$abundance_compartments) <- compartment_names
          }
        } # abundance results

        if ("ema" %in% results_selection) {
          abundance_count_min <-
            pmin(all_abundance_count, results$abundance_count_min)
          if (replicates > 1) {
            results$all$ema[tm, season] <- results$all$ema[tm, season] +
              (abundance_count_min - results$all$ema[tm, season]) / r
          } else {
            results$all$ema[tm, season] <- abundance_count_min
          }
          results$abundance_count_min <- abundance_count_min
        }

        # Extirpation results
        if ("extirpation" %in% results_selection) {
          if (replicates > 1) {
            results$extirpation[, r] <- pmin(results$extirpation[, r],
                                             rep(tm - 1 + season / seasons,
                                                 populations),
                                             na.rm = TRUE)
            results$extirpation[which(as.logical(abundance_count)), r] <- NA
          } else {
            # single replicate
            results$extirpation <- pmin(results$extirpation,
                                        rep(tm + season / seasons, populations),
                                        na.rm = TRUE)
            results$extirpation[which(as.logical(abundance_count))] <- NA
          }
          results$all$extirpation[r] <- min(results$extirpation[r],
                                            tm - 1 + season / seasons,
                                            na.rm = TRUE)
          results$all$extirpation[r][as.logical(all_abundance_count)] <- NA
        }

        # Extinction location results
        if ("extinction_location" %in% results_selection) {
          if (all_abundance_count > 0) {
            # => not extinct
            results$last_occupied_abundance_count <- abundance_count
          }
        }

        # Occupancy results
        if ("occupancy" %in% results_selection) {
          if (replicates > 1) {
            if ("replicate" %in% results_selection) {
              results$all$occupancy[tm, season, r] <-
                sum(as.logical(abundance_count))
            } else {
              # summarize (default)
              occupancy_count <- sum(as.logical(abundance_count))
              previous_occupancy_mean <-
                results$all$occupancy$mean[tm, season]
              results$all$occupancy$mean[tm, season] <-
                previous_occupancy_mean +
                (occupancy_count - previous_occupancy_mean) / r
              previous_occupancy_var_by_r <-
                results$all$occupancy$sd[tm, season]
              results$all$occupancy$sd[tm, season] <-
                previous_occupancy_var_by_r +
                (occupancy_count - previous_occupancy_mean) *
                (occupancy_count - results$all$occupancy$mean[tm, season])
              results$all$occupancy$min[tm, season] <- ifelse(
                r == 1,
                occupancy_count,
                min(results$all$occupancy$min[tm, season],
                    occupancy_count)
              )
              results$all$occupancy$max[tm, season] <- ifelse(
                r == 1,
                occupancy_count,
                max(results$all$occupancy$max[tm, season],
                    occupancy_count)
              )
            }
          } else {
            # single replicate
            results$all$occupancy[tm, season] <- sum(as.logical(abundance_count))
          }
        }
      } # results calculated from abundance

      if (!is.null(harvested)) {
        # results calculated from harvested

        # Select harvested counts for included stages and compartments
        harvested_count <- as.array(.colSums(harvested,
                                             m = stages * compartments,
                                             n = populations))
        all_harvested_count <- sum(harvested_count)

        if (results_breakdown == "segments") {
          harvested_segments_count <-
            lapply(1:(stages * compartments), function(s) {
              as.array(.colSums(harvested[s,], m = length(s), n = populations))
            })
          all_harvested_segments_count <-
            lapply(harvested_segments_count, sum)
        } else if (results_breakdown == "stages") {
          harvested_stages_count <- lapply(1:stages, function(s) {
            stage_indices <- seq(s, nrow(harvested), by = stages)
            as.array(.colSums(
              harvested[stage_indices,],
              m = length(stage_indices),
              n = populations
            ))
          })
          all_harvested_stages_count <-
            lapply(harvested_stages_count, sum)
        } else if (results_breakdown == "compartments") {
          harvested_compartments_count <- lapply(1:compartments, function(s) {
            compartment_indices <- seq(stages * s - (stages - 1), stages * s, 1)
            as.array(.colSums(
              harvested[compartment_indices,],
              m = length(compartment_indices),
              n = populations
            ))
          })
          all_harvested_compartments_count <-
            lapply(harvested_compartments_count, sum)
        }

        # Harvested results
        if ("harvested" %in% results_selection) {
          if (replicates > 1) {
            if ("replicate" %in% results_selection) {
              results$harvested[, tm, season, r] <- harvested_count
              results$all$harvested[tm, season, r] <-
                all_harvested_count
            } else {
              # summarize (default)
              previous_harvested_mean <-
                results$harvested$mean[, tm, season]
              results$harvested$mean[, tm, season] <-
                previous_harvested_mean +
                (harvested_count - previous_harvested_mean) / r
              previous_harvested_var_by_r <-
                results$harvested$sd[, tm, season]
              results$harvested$sd[, tm, season] <-
                previous_harvested_var_by_r +
                (harvested_count - previous_harvested_mean) *
                (harvested_count - results$harvested$mean[, tm, season])

              if (r == 1) {
                results$harvested$min[, tm, season] <- harvested_count
                results$harvested$max[, tm, season] <- harvested_count
              } else {
                results$harvested$min[, tm, season] <- pmin(results$harvested$min[, tm, season],
                                                            harvested_count)
                results$harvested$max[, tm, season] <- pmax(results$harvested$max[, tm, season],
                                                            harvested_count)
              }

              previous_all_harvested_mean <-
                results$all$harvested$mean[tm, season]
              results$all$harvested$mean[tm, season] <-
                previous_all_harvested_mean +
                (all_harvested_count - previous_all_harvested_mean) / r
              previous_all_harvested_var_by_r <-
                results$all$harvested$sd[tm, season]
              results$all$harvested$sd[tm, season] <-
                previous_all_harvested_var_by_r +
                (all_harvested_count - previous_all_harvested_mean) *
                (all_harvested_count - results$all$harvested$mean[tm, season])
              results$all$harvested$min[tm, season] <- ifelse(
                r > 1,
                min(results$all$harvested$min[tm, season],
                    all_harvested_count),
                all_harvested_count
              )
              results$all$harvested$max[tm, season] <- ifelse(
                r > 1,
                max(results$all$harvested$max[tm, season],
                    all_harvested_count),
                all_harvested_count
              )
            }
          } else {
            # single replicate
            results$harvested[, tm, season] <- harvested_count
            results$all$harvested[tm, season] <- all_harvested_count
          }

          if (results_breakdown == "segments") {
            segment_names <- names(results$harvested_segments)
            for (i in 1:(stages * compartments)) {
              if (replicates > 1) {
                if ("replicate" %in% results_selection) {
                  results$harvested_segments[[i]][, tm, season, r] <-
                    harvested_segments_count[[i]]
                  results$all$harvested_segments[[i]][tm, season, r] <-
                    all_harvested_segments_count[[i]]
                  names(results$all$harvested_segments) <- segment_names
                } else {
                  previous_harvested_segments_mean <-
                    results$harvested_segments[[i]]$mean[, tm, season]
                  results$harvested_segments[[i]]$mean[, tm, season] <-
                    previous_harvested_segments_mean +
                    (harvested_segments_count[[i]] - previous_harvested_segments_mean) / r
                  previous_harvested_segments_var_by_r <-
                    results$harvested_segments[[i]]$sd[, tm, season]
                  results$harvested_segments[[i]]$sd[, tm, season] <-
                    previous_harvested_segments_var_by_r +
                    (harvested_segments_count[[i]] - previous_harvested_segments_mean) *
                    (harvested_segments_count[[i]] - results$harvested_segments[[i]]$mean[, tm, season])
                  names(results$harvested_segments) <- segment_names
                  if (r == 1) {
                    results$harvested_segments[[i]]$min[, tm, season] <-
                      harvested_segments_count[[i]]
                    results$harvested_segments[[i]]$max[, tm, season] <-
                      harvested_segments_count[[i]]
                  } else {
                    results$harvested_segments[[i]]$min[, tm, season] <- pmin(results$harvested_segments[[i]]$min[, tm, season],
                                                                              harvested_segments_count[[i]])
                    results$harvested_segments[[i]]$max[, tm, season] <-
                      pmax(results$harvested_segments[[i]]$max[, tm, season],
                           harvested_segments_count[[i]])
                  }

                  previous_all_harvested_segments_mean <-
                    results$all$harvested_segments[[i]]$mean[tm, season]
                  results$all$harvested_segments[[i]]$mean[tm, season] <-
                    previous_all_harvested_segments_mean +
                    (
                      all_harvested_segments_count[[i]] - previous_all_harvested_segments_mean
                    ) / r
                  previous_all_harvested_segments_var_by_r <-
                    results$all$harvested_segments[[i]]$sd[tm, season]
                  results$all$harvested_segments[[i]]$sd[tm, season] <-
                    previous_all_harvested_segments_var_by_r +
                    (
                      all_harvested_segments_count[[i]] - previous_all_harvested_segments_mean
                    ) *
                    (
                      all_harvested_segments_count[[i]] - results$all$harvested_segments[[i]]$mean[tm, season]
                    )

                  results$all$harvested_segments[[i]]$min[tm, season] <-
                    ifelse(
                      r > 1,
                      min(
                        results$all$harvested_segments[[i]]$min[tm, season],
                        all_harvested_segments_count[[i]]
                      ),
                      all_harvested_segments_count[[i]]
                    )
                  results$all$harvested_segments[[i]]$max[tm, season] <-
                    ifelse(
                      r > 1,
                      max(
                        results$all$harvested_segments[[i]]$max[tm, season],
                        all_harvested_segments_count[[i]]
                      ),
                      all_harvested_segments_count[[i]]
                    )
                  names(results$all$harvested_segments) <- segment_names
                }
              } else {
                # single replicate
                results$harvested_segments[[i]][, tm, season] <-
                  harvested_segments_count[[i]]
                results$all$harvested_segments[[i]][tm, season] <-
                  all_harvested_segments_count[[i]]
                names(results$all$harvested_segments) <- segment_names
              }
            }
          }

          if (results_breakdown == "stages") {
            stage_names <- names(results$harvested_stages)
            for (i in 1:stages) {
              if (replicates > 1) {
                if ("replicate" %in% results_selection) {
                  results$harvested_stages[[i]][, tm, season, r] <-
                    harvested_stages_count[[i]]
                  results$all$harvested_stages[[i]][tm, season, r] <-
                    all_harvested_stages_count[[i]]
                  names(results$all$harvested_stages) <- stage_names
                  names(results$harvested_stages) <- stage_names
                } else {
                  previous_harvested_stages_mean <-
                    results$harvested_stages[[i]]$mean[, tm, season]
                  results$harvested_stages[[i]]$mean[, tm, season] <-
                    previous_harvested_stages_mean +
                    (harvested_stages_count[[i]] - previous_harvested_stages_mean) / r
                  previous_harvested_stages_var_by_r <-
                    results$harvested_stages[[i]]$sd[, tm, season]
                  results$harvested_stages[[i]]$sd[, tm, season] <-
                    previous_harvested_stages_var_by_r +
                    (harvested_stages_count[[i]] - previous_harvested_stages_mean) *
                    (harvested_stages_count[[i]] - results$harvested_stages[[i]]$mean[, tm, season])
                  names(results$harvested_stages) <- stage_names
                  if (r == 1) {
                    results$harvested_stages[[i]]$min[, tm, season] <-
                      harvested_stages_count[[i]]
                    results$harvested_stages[[i]]$max[, tm, season] <-
                      harvested_stages_count[[i]]
                  } else {
                    results$harvested_stages[[i]]$min[, tm, season] <- pmin(results$harvested_stages[[i]]$min[, tm, season],
                                                                            harvested_stages_count[[i]])
                    results$harvested_stages[[i]]$max[, tm, season] <-
                      pmax(results$harvested_stages[[i]]$max[, tm, season],
                           harvested_stages_count[[i]])
                  }

                  previous_all_harvested_stages_mean <-
                    results$all$harvested_stages[[i]]$mean[tm, season]
                  results$all$harvested_stages[[i]]$mean[tm, season] <-
                    previous_all_harvested_stages_mean +
                    (all_harvested_stages_count[[i]] - previous_all_harvested_stages_mean) / r
                  previous_all_harvested_stages_var_by_r <-
                    results$all$harvested_stages[[i]]$sd[tm, season]
                  results$all$harvested_stages[[i]]$sd[tm, season] <-
                    previous_all_harvested_stages_var_by_r +
                    (all_harvested_stages_count[[i]] - previous_all_harvested_stages_mean) *
                    (
                      all_harvested_stages_count[[i]] - results$all$harvested_stages[[i]]$mean[tm, season]
                    )

                  results$all$harvested_stages[[i]]$min[tm, season] <-
                    ifelse(
                      r > 1,
                      min(
                        results$all$harvested_stages[[i]]$min[tm, season],
                        all_harvested_stages_count[[i]]
                      ),
                      all_harvested_stages_count[[i]]
                    )
                  results$all$harvested_stages[[i]]$max[tm, season] <-
                    ifelse(
                      r > 1,
                      max(
                        results$all$harvested_stages[[i]]$max[tm, season],
                        all_harvested_stages_count[[i]]
                      ),
                      all_harvested_stages_count[[i]]
                    )
                  names(results$all$harvested_stages) <- stage_names
                }
              } else {
                # single replicate
                results$harvested_stages[[i]][, tm, season] <-
                  harvested_stages_count[[i]]
                results$all$harvested_stages[[i]][tm, season] <-
                  all_harvested_stages_count[[i]]
                names(results$all$harvested_stages) <- stage_names
              }
            }
          }

          if (results_breakdown == "compartments") {
            compartment_names <- names(results$harvested_compartments)
            for (i in 1:compartments) {
              if (replicates > 1) {
                if ("replicate" %in% results_selection) {
                  results$harvested_compartments[[i]][, tm, season, r] <-
                    harvested_compartments_count[[i]]
                  results$all$harvested_compartments[[i]][tm, season, r] <-
                    all_harvested_compartments_count[[i]]
                  names(results$all$harvested_compartments) <- compartment_names
                  names(results$harvested_compartments) <- compartment_names
                } else {
                  previous_harvested_compartments_mean <-
                    results$harvested_compartments[[i]]$mean[, tm, season]
                  results$harvested_compartments[[i]]$mean[, tm, season] <-
                    previous_harvested_compartments_mean +
                    (
                      harvested_compartments_count[[i]] - previous_harvested_compartments_mean
                    ) / r
                  previous_harvested_compartments_var_by_r <-
                    results$harvested_compartments[[i]]$sd[, tm, season]
                  results$harvested_compartments[[i]]$sd[, tm, season] <-
                    previous_harvested_compartments_var_by_r +
                    (
                      harvested_compartments_count[[i]] - previous_harvested_compartments_mean
                    ) *
                    (
                      harvested_compartments_count[[i]] - results$harvested_compartments[[i]]$mean[, tm, season]
                    )
                  names(results$harvested_compartments) <- compartment_names

                  if (r == 1) {
                    results$harvested_compartments[[i]]$min[, tm, season] <-
                      harvested_compartments_count[[i]]
                    results$harvested_compartments[[i]]$max[, tm, season] <-
                      harvested_compartments_count[[i]]
                  } else {
                    results$harvested_compartments[[i]]$min[, tm, season] <- pmin(
                      results$harvested_compartments[[i]]$min[, tm, season],
                      harvested_compartments_count[[i]]
                    )
                    results$harvested_compartments[[i]]$max[, tm, season] <-
                      pmax(
                        results$harvested_compartments[[i]]$max[, tm, season],
                        harvested_compartments_count[[i]]
                      )
                  }

                  previous_all_harvested_compartments_mean <-
                    results$all$harvested_compartments[[i]]$mean[tm, season]
                  results$all$harvested_compartments[[i]]$mean[tm, season] <-
                    previous_all_harvested_compartments_mean +
                    (
                      all_harvested_compartments_count[[i]] - previous_all_harvested_compartments_mean
                    ) / r
                  previous_all_harvested_compartments_var_by_r <-
                    results$all$harvested_compartments[[i]]$sd[tm, season]
                  results$all$harvested_compartments[[i]]$sd[tm, season] <-
                    previous_all_harvested_compartments_var_by_r +
                    (
                      all_harvested_compartments_count[[i]] - previous_all_harvested_compartments_mean
                    ) *
                    (
                      all_harvested_compartments_count[[i]] - results$all$harvested_compartments[[i]]$mean[tm, season]
                    )

                  results$all$harvested_compartments[[i]]$min[tm, season] <-
                    ifelse(
                      r > 1,
                      min(
                        results$all$harvested_compartments[[i]]$min[tm, season],
                        all_harvested_compartments_count[[i]]
                      ),
                      all_harvested_compartments_count[[i]]
                    )
                  results$all$harvested_compartments[[i]]$max[tm, season] <-
                    ifelse(
                      r > 1,
                      max(
                        results$all$harvested_compartments[[i]]$max[tm, season],
                        all_harvested_compartments_count[[i]]
                      ),
                      all_harvested_compartments_count[[i]]
                    )
                  names(results$all$harvested_compartments) <- compartment_names
                }
              } else {
                # single replicate
                results$harvested_compartments[[i]][, tm, season] <-
                  harvested_compartments_count[[i]]
                results$all$harvested_compartments[[i]][tm, season] <-
                  all_harvested_compartments_count[[i]]
                names(results$all$harvested_compartments) <- compartment_names
              }
            }
          }
        }
      } # results calculated from harvested

      return(results)
    } # End calculate at season

  ## Create a nested function for calculating/collecting (non-NULL) results at replicate ##
  calculate_at_replicate <-
    function(r, segment_abundance, results) {
      # Sum abundance counts for included segments
      all_abundance_count <- sum(.colSums(
        segment_abundance,
        m = stages * compartments,
        n = populations
      ))

      # Extinction location results
      if ("extinction_location" %in% results_selection &&
          !is.null(coordinates)) {
        if (all_abundance_count == 0) {
          # => extinct
          last_pop_indices <-
            which(as.logical(results$last_occupied_abundance_count))
          if (length(last_pop_indices) > 1) {
            # abundance-weighted average of last occupied coordinates
            abundance_weights <-
              matrix(rep(results$last_occupied_abundance_count[last_pop_indices], 2),
                     ncol = 2)
            results$all$extinction_location[r,] <-
              (
                .colSums(
                  coordinates[last_pop_indices,] * abundance_weights,
                  m = length(last_pop_indices),
                  n = 2
                ) / .colSums(
                  abundance_weights,
                  m = length(last_pop_indices),
                  n = 2
                )
              )
          } else {
            # last occupied coordinates
            results$all$extinction_location[r,] <-
              as.numeric(coordinates[last_pop_indices,])
          }
        }
      }

      return(results)
    }

  # Create a nested function for finalizing result calculations at the end of all replicates
  finalize_attributes <- function(results) {
    # Finalize/summarize results
    if (replicates > 1 && !("replicate" %in% results_selection)) {
      if ("abundance" %in% results_selection) {
        results$abundance$sd <-
          sqrt(results$abundance$sd / (replicates - 1))
        results$all$abundance$sd <-
          sqrt(results$all$abundance$sd / (replicates - 1))

        if (results_breakdown == "segments") {
          results$abundance_segments <-
            lapply(results$abundance_segments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
          results$all$abundance_segments <-
            lapply(results$all$abundance_segments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
        } else if (results_breakdown == "stages") {
          results$abundance_stages <-
            lapply(results$abundance_stages, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
          results$all$abundance_stages <-
            lapply(results$all$abundance_stages, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
        } else if (results_breakdown == "compartments") {
          results$abundance_compartments <-
            lapply(results$abundance_compartments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
          results$all$abundance_compartments <-
            lapply(results$all$abundance_compartments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
        }
      }

      if ("harvested" %in% results_selection) {
        results$harvested$sd <-
          sqrt(results$harvested$sd / (replicates - 1))
        results$all$harvested$sd <-
          sqrt(results$all$harvested$sd / (replicates - 1))

        if (results_breakdown == "segments") {
          results$harvested_segments <-
            lapply(results$harvested_segments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
          results$all$harvested_segments <-
            lapply(results$all$harvested_segments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
        } else if (results_breakdown == "stages") {
          results$harvested_stages <-
            lapply(results$harvested_stages, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
          results$all$harvested_stages <-
            lapply(results$all$harvested_stages, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
        } else if (results_breakdown == "compartments") {
          results$harvested_compartments <-
            lapply(results$harvested_compartments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
          results$all$harvested_compartments <-
            lapply(results$all$harvested_compartments, function(x) {
              x$sd <- sqrt(x$sd / (replicates - 1))
              return(x)
            })
        }
      }

      if ("occupancy" %in% results_selection) {
        results$all$occupancy$sd <-
          sqrt(results$all$occupancy$sd / (replicates - 1))
      }

      # Summarize extirpation dependent on the presence of NAs
      # (since mean and sd cannot be computed with NAs => no extirpation)
      if ("extirpation" %in% results_selection) {
        if (any(is.na(results$extirpation))) {
          # 5 number summary when NAs present
          results$extirpation[which(is.na(results$extirpation))] <-
            time_steps + 1
          results$extirpation <-
            apply(results$extirpation, 1, stats::fivenum)
          results$extirpation[which(results$extirpation > time_steps)] <-
            NA
          results$extirpation <- list(
            min = results$extirpation[1,],
            q1 = results$extirpation[2,],
            median = results$extirpation[3,],
            q3 = results$extirpation[4,],
            max = results$extirpation[5,]
          )
        } else {
          # mean, sd, min, max
          results$extirpation <- list(
            mean = apply(results$extirpation, 1, mean),
            sd = apply(results$extirpation, 1, stats::sd),
            min = apply(results$extirpation, 1, min),
            max = apply(results$extirpation, 1, max)
          )
        }
      }
    }

    # Remove working abundance counts
    results$abundance_count_min <- NULL
    results$last_occupied_abundance_count <- NULL

    return(results)
  } # End finalize attributes

  return(
    list(
      initialize_attributes = initialize_attributes,
      initialize_replicate = initialize_replicate,
      calculate_at_season = calculate_at_season,
      calculate_at_replicate = calculate_at_replicate,
      finalize_attributes = finalize_attributes
    )
  )
} # End disease results function
