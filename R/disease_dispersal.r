#' Nested functions for population dispersal in a disease ecology simulation.
#'
#' Modular functions for the disease simulator for performing dispersal of
#' segment (stage by compartment) abundance at a specified time step via
#' dispersal rates provided. Dispersal can be handled identically for all stages
#' and compartments, or can be handled differently by stage, compartment, or
#' both.
#'
#' @param replicates Number of replicate simulation runs.
#' @param time_steps Number of simulation time steps.
#' @param populations Number of populations.
#' @param demographic_stochasticity Boolean for optionally choosing demographic
#'   stochasticity for the transformation.
#' @param dispersal Must be a list. Either a list of matrices of dispersal rates
#'   between populations (source columns to target rows) or a list of data
#'   frames of non-zero dispersal rates and indices for constructing a compact
#'   dispersal matrix, and optionally a list of lists of data frames showing
#'   changing rates over time (as per class [`poems::DispersalGenerator`]
#'   \emph{dispersal_data} attribute). Alternatively a list of user-defined
#'   functions may be used: \code{function(params)}, where \emph{params} is a
#'   list passed to the function containing:
#'   \describe{
#'     \item{\code{replicates}}{Number of replicate simulation runs.}
#'     \item{\code{time_steps}}{Number of simulation time steps.}
#'     \item{\code{populations}}{Number of populations.}
#'     \item{\code{stages}}{Number of life cycle stages.}
#'     \item{\code{compartments}}{Number of disease compartments.}
#'     \item{\code{dispersal_type}}{Must be "pooled", "stages", "compartments", or "segments". This indicates whether dispersal should be handled differently by stage and/or compartment.}
#'     \item{\code{demographic_stochasticity}}{Boolean for optionally choosing demographic stochasticity for the transformation.}
#'     \item{\code{dispersal_source_n_k}}{Dispersal proportion (p) density dependence via source population abundance divided by carrying capacity (n/k), where p is reduced via a linear slope (defined by two list items) from n/k <= \emph{cutoff} (p = 0) to n/k >= \emph{threshold}.}
#'     \item{\code{dispersal_target_k}}{Dispersal rate (r) density dependence via target population carrying capacity (k), where r is reduced via a linear slope (through the origin) when k <= \emph{threshold}.}
#'     \item{\code{dispersal_target_n}}{Dispersal rate (r) density dependence via target population abundance (n), where r is reduced via a linear slope (defined by two list items) from n >= \emph{threshold} to n <= \emph{cutoff} (r = 0) or vice versa.}
#'     \item{\code{dispersal_target_n_k}}{Dispersal rate (r) density dependence via target population abundance divided by carrying capacity (n/k), where r is reduced via a linear slope (defined by two list items) from n/k >= \emph{threshold} to n/k <= \emph{cutoff} (r = 0) or vice versa.}
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for each population at time step.}
#'     \item{\code{segment_abundance}}{Matrix of abundance for each stage by compartment (rows) and population (columns) at time step.}
#'     \item{\code{simulator}}{[`poems::SimulatorReference`] object with dynamically accessible \emph{attached} and \emph{results} lists.}
#'   }
#'   returns the post-dispersal abundance matrix
#' @param dispersal_type Must be "pooled", "stages", "compartments", or
#'   "segments". This indicates whether dispersal should be handled differently
#'   by stage and/or compartment.
#' @param dispersal_source_n_k Dispersal proportion (p) density dependence via
#'   source population abundance divided by carrying capacity (n/k), where p is
#'   reduced via a linear slope (defined by two list items) from n/k <=
#'   \emph{cutoff} (p = 0) to n/k >= \emph{threshold} or vice versa.
#' @param dispersal_target_k Dispersal rate (r) density dependence via target
#'   population carrying capacity (k), where r is reduced via a linear slope
#'   (through the origin) when k <= \emph{threshold}.
#' @param dispersal_target_n Dispersal rate (r) density dependence via target
#'   population abundance (n), where r is reduced via a linear slope (defined by
#'   two list items) from n >= \emph{threshold} to n <= \emph{cutoff} (r = 0) or
#'   visa-versa.
#' @param dispersal_target_n_k Dispersal rate (r) density dependence via target
#'   population abundance divided by carrying capacity (n/k), where r is reduced
#'   via a linear slope (defined by two list items) from n/k >= \emph{threshold}
#'   to n/k <= \emph{cutoff} (r = 0) or vice versa.
#' @param stages Number of life cycle stages.
#' @param compartments Number of disease compartments.
#' @param simulator [`poems::SimulatorReference`] object with dynamically
#'   accessible \emph{attached} and \emph{results} lists.
#' @return Dispersal function: \code{function(r, tm, carrying_capacity,
#'   segment_abundance)}, where:
#'   \describe{
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for each population at time step.}
#'     \item{\code{segment_abundance}}{Matrix of abundance for each stage by compartment (rows) and population (columns) at time step.}
#'     \item{\code{returns}}{New stage abundance matrix with dispersal applied.}
#'   }
#' @import purrr
#' @export disease_dispersal

disease_dispersal <- function(replicates,
                              time_steps,
                              populations,
                              demographic_stochasticity,
                              dispersal,
                              dispersal_type,
                              dispersal_source_n_k = NULL,
                              dispersal_target_k = NULL,
                              dispersal_target_n = NULL,
                              dispersal_target_n_k = NULL,
                              stages = NULL,
                              compartments = NULL,
                              simulator) {
  if (is.null(dispersal)) { # no dispersal
    return(NULL)
  }

  # Sanity checks
  if (dispersal_type == "pooled" && length(dispersal) != 1) {
      cli_abort(c("Error: Dispersal length should be 1 for 'pooled' dispersal
                  type.",
                "x" = "Dispersal length is {length(dispersal)}."))
  } else if (dispersal_type == "stages" && length(dispersal) != stages) {
    cli_abort(c("Error: Dispersal length should be equal to the number of
    stages for 'stages' dispersal type.",
              "x" = "There are {stages} stages and dispersal length is
              {length(dispersal)}."))
  } else if (dispersal_type == "compartments" &&
              length(dispersal) != compartments) {
    cli_abort(c("Error: Dispersal length should be equal to the number of
    compartments for 'compartments' dispersal type.",
              "x" = "There are {compartments} compartments and dispersal
              length is {length(dispersal)}."))
  } else if (dispersal_type == "segments" &&
              length(dispersal) != stages * compartments) {
    cli_abort(c("Error: Dispersal length should be equal to the number of
    stages multiplied by the number of compartments for 'segments' dispersal
                type.",
    "x" = "There are {stages} stages and {compartments} compartments and
    dispersal length is {length(dispersal)}."))
  }

  # User-defined function?
  if (length(which(unlist(lapply(dispersal, is.function))))) {
    # List of parameters to pass to the user-defined function
    params <- list(replicates = replicates, time_steps = time_steps,
                   populations = populations, stages = stages,
                   compartments = compartments, dispersal_type = dispersal_type,
                   demographic_stochasticity = demographic_stochasticity,
                   dispersal_source_n_k = dispersal_source_n_k,
                   dispersal_target_k = dispersal_target_k,
                   dispersal_target_n = dispersal_target_n,
                   dispersal_target_n_k = dispersal_target_n_k,
                   simulator = simulator)

    ## Create a nested function for applying user-defined dispersal of segment abundance ##
    dispersal_function <- function(r, tm, carrying_capacity, segment_abundance) {

      if (dispersal_type == "stages") {
        step_indices <- lapply(1:stages, function(s) {
          step_indices <- seq(s, nrow(segment_abundance), by = stages)
        })
      } else if (dispersal_type == "compartments") {
        step_indices <- lapply(1:compartments, function(s) {
          step_indices <- seq(stages * s - (stages - 1), stages * s, 1)
        })
      } else if (dispersal_type == "segments") {
        step_indices_vector <- c(1:(stages*compartments))
        step_indices <- as.list(step_indices_vector)
      } else if (dispersal_type == "pooled") {
        step_indices <- rep(list(1:(stages*compartments)))
      }

      for (f in 1:length(dispersal)) {

        # Add attributes to be made available to the user-defined function
        params$r <- r
        params$tm <- tm
        params$carrying_capacity <- carrying_capacity
        params$segment_abundance <- segment_abundance[step_indices[[f]], , drop = F]
        params$occupied_indices <- which(apply(params$segment_abundance, 2, sum) > 0)

        # Run user-defined dispersal function
        tryCatch({
          segment_abundance[step_indices[[f]],] <- dispersal[[f]](params)
        },
        error = function(e){
          stop(paste("Error produced within user-defined dispersal function:", as.character(e)), call. = FALSE)
        })

        # Warn if any negative or non-finite
        if (any(!is.finite(segment_abundance))) {
          warning("Non-finite abundances returned by user-defined dispersal function", call. = FALSE)
        }
        if (any(segment_abundance[which(is.finite(segment_abundance))] < 0)) {
          warning("Negative abundances returned by user-defined dispersal function", call. = FALSE)
        }
      }

      return(segment_abundance)
    }
    return(dispersal_function)
  } # end if function list

  # Initialize reusable dispersal attributes
  if (length(which(unlist(lapply(dispersal, is.matrix))))) {

    # Initialize lists
    dispersal_data_list <- list()
    dispersal_compact_rows_list <- list()
    dispersals_change_over_time_list <- list()
    dispersal_stages <- map_lgl(dispersal, ~nrow(.x) > 0)

    # Loop over each matrix in the dispersal list
    for (i in seq_along(dispersal[dispersal_stages])) {
      if (is.matrix(dispersal[[i]])) { # create compact matrix data

        # Calculate the indices of non-zero dispersals
        dispersal_data <- data.frame(which(dispersal[[i]] > 0, arr.ind = TRUE))
        dispersal_data <- dispersal_data[order(dispersal_data[, 2], dispersal_data[, 1]),]
        names(dispersal_data) <- c("target_pop", "source_pop")

        # Calculate indices for constructing compacted dispersal matrices for emigrants and immigrants
        dispersal_rows <- tabulate(dispersal_data$source_pop, nbins = populations)
        dispersal_cols <- tabulate(dispersal_data$target_pop, nbins = populations)
        dispersal_compact_rows <- max(dispersal_rows, dispersal_cols)
        compact_emigrant_matrix <- array(1:dispersal_compact_rows, c(dispersal_compact_rows, populations))
        compact_immigrant_matrix <- compact_emigrant_matrix*(compact_emigrant_matrix <= matrix(dispersal_cols, nrow = dispersal_compact_rows, ncol = populations, byrow = TRUE))
        compact_emigrant_matrix <- compact_emigrant_matrix*(compact_emigrant_matrix <= matrix(dispersal_rows, nrow = dispersal_compact_rows, ncol = populations, byrow = TRUE))

        # Map the row of each compact matrix to the original target (for emigrants) or source (for immigrants) populations
        dispersal_data$emigrant_row <- which(compact_emigrant_matrix > 0, arr.ind = TRUE, useNames = FALSE)[,1]
        dispersal_data$immigrant_row <- which(compact_immigrant_matrix > 0, arr.ind = TRUE, useNames = FALSE)[,1]
        target_sorted_indices <- dispersal_data[order(dispersal_data$target_pop, dispersal_data$source_pop), c("target_pop", "source_pop")]
        dispersal_data$immigrant_row <- dispersal_data$immigrant_row[order(target_sorted_indices$source_pop, target_sorted_indices$target_pop)]

        # Add dispersal rates
        dispersal_data$dispersal_rate <- dispersal[[i]][as.matrix(dispersal_data[c("target_pop", "source_pop")])]
        dispersals_change_over_time <- FALSE

        # Add to lists
        dispersal_data_list[[i]] <- dispersal_data
        dispersal_compact_rows_list[[i]] <- dispersal_compact_rows
        dispersals_change_over_time_list[[i]] <- dispersals_change_over_time

        # Release variables from memory
        dispersal_rows <- NULL; dispersal_cols <- NULL; compact_emigrant_matrix = NULL; compact_immigrant_matrix = NULL; target_sorted_indices = NULL
      }
    }
  } else if (is.list(dispersal) && all(sapply(dispersal, is.list)) &&
  all(sapply(unlist(dispersal, recursive = FALSE), is.data.frame)) &&
  nrow(dispersal[[1]][[1]]) > 0) {
    # Initialize lists
    dispersal_data_list <- list()
    dispersal_compact_rows_list <- list()
    dispersals_change_over_time_list <- list()
    dispersal_data_changes_list <- list()

    if (!is.list(dispersal[[1]])) {
      dispersal <- list(dispersal)
    }

    dispersal_stages <- !map_lgl(dispersal, ~all(map_lgl(.x, ~nrow(.x) == 0)))

    # Loop over each list in the dispersal list
    for (i in seq_along(dispersal[dispersal_stages])) {
      # Unpack dispersal data and determine compact matrix dimensions
      dispersal_data <- dispersal[[i]][[1]]
      dispersal_compact_rows <- max(dispersal_data[, c("emigrant_row", "immigrant_row")])

      # Are dispersals changing over time?
      dispersals_change_over_time <- (length(dispersal[[i]]) > 1)
      if (dispersals_change_over_time) {
        dispersal_data_changes <- dispersal[[i]]
        dispersal_data_changes[[1]] <- dispersal_data_changes[[1]][NULL,]
      }

      # Add to lists
      dispersal_data_list[[i]] <- dispersal_data
      dispersal_compact_rows_list[[i]] <- dispersal_compact_rows
      dispersals_change_over_time_list[[i]] <- dispersals_change_over_time
      dispersal_data_changes_list[[i]] <- if (exists("dispersal_data_changes")) dispersal_data_changes else NULL
    }
  } else {
    return(NULL)
  }

  # Release dispersal from memory
  dispersal <- NULL

  # Create a function to generate a compact matrix for each element in the lists
  generate_compact_matrix <- function(dispersal_data, dispersal_compact_rows) {
    dispersal_compact_matrix <- array(0, c(dispersal_compact_rows, populations))
    dispersal_compact_matrix[as.matrix(dispersal_data[, c("emigrant_row", "source_pop")])] <- dispersal_data$dispersal_rate
    return(dispersal_compact_matrix)
  }

  # Use lapply to apply the function to each element in the lists
  dispersal_compact_matrix_list <- mapply(generate_compact_matrix, dispersal_data_list,
                                          dispersal_compact_rows_list, SIMPLIFY=FALSE)

  # Does dispersal depend on source population abundance N divided by carrying capacity K?
  dispersal_depends_on_source_pop_n_k <- (is.list(dispersal_source_n_k) && (is.numeric(dispersal_source_n_k$cutoff) ||
                                                                              is.numeric(dispersal_source_n_k$threshold)))

  # Does dispersal depend on target population carrying capacity K, abundance N, or N/K?
  dispersal_depends_on_target_pop_k <- is.numeric(dispersal_target_k)
  dispersal_depends_on_target_pop_n <- (is.list(dispersal_target_n) && (is.numeric(dispersal_target_n$threshold) ||
                                                                          is.numeric(dispersal_target_n$cutoff)))
  dispersal_depends_on_target_pop_n_k <- (is.list(dispersal_target_n_k) && (is.numeric(dispersal_target_n_k$threshold) ||
                                                                              is.numeric(dispersal_target_n_k$cutoff)))

  # Setup density dependence dispersal parameters
  if (dispersal_depends_on_source_pop_n_k) {

    # Convert NULL to zero in source N/K cutoff or one in threshold
    if (dispersal_depends_on_source_pop_n_k) {
      if (is.null(dispersal_source_n_k$cutoff)) dispersal_source_n_k$cutoff <- 0
      if (is.null(dispersal_source_n_k$threshold)) dispersal_source_n_k$threshold <- 1
    }

    # Check threshold > cutoff
    if (dispersal_source_n_k$threshold <= dispersal_source_n_k$cutoff) {
      dispersal_depends_on_source_pop_n_k <- FALSE
      cli_warn("Dispersal density dependence for source N/K threshold must be greater than cutoff.",
               "i" = "Source threshold is {dispersal_source_n_k$threshold} and source cutoff is
                     {dispersal_source_n_k$cutoff}.",
               "x" = "Dispersal density dependence for source N/K not used.")
    }
  }

  if (dispersal_depends_on_target_pop_k || dispersal_depends_on_target_pop_n || dispersal_depends_on_target_pop_n_k) {

    if (dispersal_depends_on_target_pop_n) {

      # Convert NULL to zero in target N threshold or cutoff
      if (is.null(dispersal_target_n$threshold)) dispersal_target_n$threshold <- 0
      if (is.null(dispersal_target_n$cutoff)) dispersal_target_n$cutoff <- 0
    }
    if (dispersal_depends_on_target_pop_n_k) {

      # Convert NULL to zero in target N/K threshold or cutoff
      if (is.null(dispersal_target_n_k$threshold)) dispersal_target_n_k$threshold <- 0
      if (is.null(dispersal_target_n_k$cutoff)) dispersal_target_n_k$cutoff <- 0
    }

    # Create a map of compact array indices for mapping dispersers (emigrants) to target populations
    # Create a function to generate a target population map for each element in the lists
    generate_target_pop_map <- function(dispersal_data, dispersal_compact_rows) {
      dispersal_target_pop_map <- array(0, c(dispersal_compact_rows, populations))
      dispersal_target_pop_map[as.matrix(dispersal_data[, c("emigrant_row", "source_pop")])] <- dispersal_data$target_pop
      return(dispersal_target_pop_map)
    }

    # Use mapply to apply the function to each element in the lists
    dispersal_target_pop_map_list <- mapply(generate_target_pop_map, dispersal_data_list,
                                            dispersal_compact_rows_list, SIMPLIFY=FALSE)
  }

  # Create a function to generate a map of compact array indices for each element in the lists
  generate_immigrant_map <- function(dispersal_data, dispersal_compact_rows) {
    dispersal_compact_indices <- array(1:(dispersal_compact_rows*populations), c(dispersal_compact_rows, populations))
    dispersal_immigrant_map <- array(0, c(dispersal_compact_rows, populations))
    dispersal_immigrant_map[as.matrix(dispersal_data[, c("emigrant_row", "source_pop")])] <- dispersal_compact_indices[as.matrix(dispersal_data[, c("immigrant_row", "target_pop")])]
    return(dispersal_immigrant_map)
  }

  # Use mapply to apply the function to each element in the lists
  dispersal_immigrant_map_list <- mapply(generate_immigrant_map, dispersal_data_list, dispersal_compact_rows_list, SIMPLIFY=FALSE)

  # Release variables from memory
  dispersal_data_list <- NULL; dispersal_compact_indices_list <- NULL

  ## Create a nested function for performing dispersal ##
  dispersal_function = function(r, tm, carrying_capacity, segment_abundance) {

    if (dispersal_type == "stages") {
      step_indices <- lapply(1:stages, function(s) {
        step_indices <- seq(s, nrow(segment_abundance), by = stages)
      })
    } else if (dispersal_type == "compartments") {
      step_indices <- lapply(1:compartments, function(s) {
        step_indices <- seq(stages * s - (stages - 1), stages * s, 1)
      })
    } else if (dispersal_type == "segments") {
      step_indices_vector <- c(1:(stages*compartments))
      step_indices <- as.list(step_indices_vector)
    } else if (dispersal_type == "pooled") {
      step_indices <- rep(list(1:(stages*compartments)))
    }

    # Calculate occupied indices
    occupied_indices_list <- map(1:nrow(segment_abundance), \(r) {
      which(apply(segment_abundance[r, , drop = F], 2, sum) > 0)
    })

    # Calculate occupied population number
    occupied_population_list <- map(occupied_indices_list, length)

    expand_lists <- function(lists, step_indices) {
      expanded_lists <- vector("list", length(unlist(step_indices)))
      for (i in seq_along(lists)) {
        expanded_lists[step_indices[[i]]] <- list(lists[[i]])
      }
      return(expanded_lists)
    }

    dispersal_compact_matrix_list <- expand_lists(dispersal_compact_matrix_list, step_indices)
    dispersals_change_over_time_list <- expand_lists(dispersals_change_over_time_list, step_indices)
    if (exists("dispersal_data_changes_list")) {
      dispersal_data_changes_list <- expand_lists(dispersal_data_changes_list, step_indices = step_indices)
    }
    dispersal_compact_rows_list <- expand_lists(dispersal_compact_rows_list, step_indices)
    dispersal_immigrant_map_list <- expand_lists(dispersal_immigrant_map_list, step_indices = step_indices)
    if (exists("dispersal_target_pop_map_list")) {
      dispersal_target_pop_map_list <- expand_lists(dispersal_target_pop_map_list, step_indices = step_indices)
    }
    dispersal_stages_expanded <- expand_lists(as.list(dispersal_stages), step_indices) |> unlist()
    if (!all(map_int(occupied_indices_list, length))) {
      dispersal_stages_expanded[map_int(occupied_indices_list, length) == 0] <- FALSE
    }
    if (all(!dispersal_stages_expanded)) {
      cli_warn("No occupied populations capable of dispersing at timestep {tm}.",
               "i" = "Dispersal not applied.")
      return(segment_abundance)
    }

    dispersal_compact_matrix_tm_list <- simulator$attached$dispersal_compact_matrix_tm_list

    # Apply any spatio-temporal dispersal changes
    apply_dispersal_changes <- function(dispersal_compact_matrix,
                                        dispersals_change_over_time,
                                        dispersal_data_changes,
                                        dispersal_compact_matrix_tm,
                                        dispersal_stages,
                                        tm) {
      if (dispersal_stages) {
        if (tm == 1 || !dispersals_change_over_time) {
          dispersal_compact_matrix_tm <- dispersal_compact_matrix
        } else if (dispersals_change_over_time &&
                   nrow(dispersal_data_changes[[tm]]) &&
                   !is.null(dispersal_compact_matrix_tm)) {
          # and tm > 1
          dispersal_compact_matrix_tm[as.matrix(dispersal_data_changes[[tm]][, c("emigrant_row", "source_pop")])] <-
            dispersal_data_changes[[tm]]$dispersal_rate
        }
      }
      return(dispersal_compact_matrix_tm)
    }

    n <- length(dispersal_compact_matrix_list)

    # Check if 'dispersal_data_changes_list' exists
    if (!exists("dispersal_data_changes_list")) {
      dispersal_data_changes_list <- vector("list", n)
    }
    if (is.null(dispersal_compact_matrix_tm_list)) {
      dispersal_compact_matrix_tm_list <- vector("list", n)
    }

    dispersal_compact_matrix_tm_list <- mapply(apply_dispersal_changes,
                                              dispersal_compact_matrix_list,
                                              dispersals_change_over_time_list,
                                              dispersal_data_changes_list,
                                              dispersal_compact_matrix_tm_list,
                                              dispersal_stages_expanded,
                                              replicate(n, list(tm)),
                                              SIMPLIFY = FALSE)

    simulator$attached$dispersal_compact_matrix_tm_list <- dispersal_compact_matrix_tm_list

    # Select dispersals for occupied populations
    occupied_dispersals_list <- dispersal_compact_matrix_tm_list |>
                                map2(occupied_indices_list, \(x, y) x[, y])

    # Calculate density abundance
    if (dispersal_depends_on_source_pop_n_k || dispersal_depends_on_target_pop_n || dispersal_depends_on_target_pop_n_k) {
      density_abundance <- colSums(segment_abundance)
    }

    # Modify dispersal rates when dispersal depends on source population N/K
    if (dispersal_depends_on_source_pop_n_k) {

      # Density dependent multipliers
      dd_multipliers <- array(1, populations)

      # Calculate the source N/K multipliers
      abundance_on_capacity <- density_abundance/carrying_capacity
      dd_multipliers[which(abundance_on_capacity <= dispersal_source_n_k$cutoff)] <- 0
      modify_pop_indices <- which(carrying_capacity > 0 & dd_multipliers > 0 &
                                    abundance_on_capacity < dispersal_source_n_k$threshold)
      dd_multipliers[modify_pop_indices] <- ((abundance_on_capacity[modify_pop_indices] -
                                                array(dispersal_source_n_k$cutoff, populations)[modify_pop_indices])/
                                               array(dispersal_source_n_k$threshold - dispersal_source_n_k$cutoff,
                                                     populations)[modify_pop_indices]*
                                               dd_multipliers[modify_pop_indices])

      # Apply modifying multipliers to dispersals
      occupied_dispersals_list <- pmap(
        list(
          occupied_dispersals_list,
          occupied_indices_list,
          dispersal_compact_rows_list,
          dispersal_stages_expanded
        ),
        \(d, i, r, l) if (l) {
          d * matrix(
            dd_multipliers[i],
            nrow = r,
            ncol = length(i),
            byrow = TRUE
          )
        } else {
          0
        }
      )

    } # dispersal depends on source pop N/K?

    # Select occupied dispersal non-zero indices
    occupied_dispersal_indices_list <- occupied_dispersals_list |>
                                       map(\(d) which(as.logical(d))) # > 0

    # Modify dispersal rates when dispersal depends on target population K, N, or N/K
    if (dispersal_depends_on_target_pop_k || dispersal_depends_on_target_pop_n || dispersal_depends_on_target_pop_n_k) {

      # Density dependent multipliers
      dd_multipliers <- array(1, populations)

      # Calculate the (below-threshold) target K multipliers
      if (dispersal_depends_on_target_pop_k) {
        modify_pop_indices <- which(carrying_capacity < dispersal_target_k)
        dd_multipliers[modify_pop_indices] <- (carrying_capacity[modify_pop_indices]/
                                                 array(dispersal_target_k, populations)[modify_pop_indices])
      }

      # Calculate the target N multipliers
      if (dispersal_depends_on_target_pop_n) {
        if (all(dispersal_target_n$threshold < dispersal_target_n$cutoff)) { # overcrowded cell avoidance \
          dd_multipliers[which(density_abundance >= dispersal_target_n$cutoff)] <- 0
          modify_pop_indices <- which(density_abundance > dispersal_target_n$threshold & dd_multipliers > 0)
          dd_multipliers[modify_pop_indices] <- ((array(dispersal_target_n$cutoff, populations)[modify_pop_indices] -
                                                    density_abundance[modify_pop_indices])/
                                                   array(dispersal_target_n$cutoff - dispersal_target_n$threshold,
                                                         populations)[modify_pop_indices]*
                                                   dd_multipliers[modify_pop_indices])
        } else if (all(dispersal_target_n$threshold > dispersal_target_n$cutoff)) { # seek company /
          dd_multipliers[which(density_abundance <= dispersal_target_n$cutoff)] <- 0
          modify_pop_indices <- which(density_abundance < dispersal_target_n$threshold & dd_multipliers > 0)
          dd_multipliers[modify_pop_indices] <- ((density_abundance[modify_pop_indices] -
                                                    array(dispersal_target_n$cutoff, populations)[modify_pop_indices])/
                                                   array(dispersal_target_n$threshold - dispersal_target_n$cutoff,
                                                         populations)[modify_pop_indices]*
                                                   dd_multipliers[modify_pop_indices])
        }
      }

      # Calculate the target N/K multipliers
      if (dispersal_depends_on_target_pop_n_k) {
        dd_multipliers[which(carrying_capacity <= 0)] <- 0
        abundance_on_capacity <- density_abundance/carrying_capacity
        if (all(dispersal_target_n_k$threshold < dispersal_target_n_k$cutoff)) { # overcrowded cell avoidance \
          dd_multipliers[which(abundance_on_capacity >= dispersal_target_n_k$cutoff)] <- 0
          modify_pop_indices <- which(abundance_on_capacity > dispersal_target_n_k$threshold & dd_multipliers > 0)
          dd_multipliers[modify_pop_indices] <- ((array(dispersal_target_n_k$cutoff, populations)[modify_pop_indices] -
                                                    abundance_on_capacity[modify_pop_indices])/
                                                   array(dispersal_target_n_k$cutoff - dispersal_target_n_k$threshold,
                                                         populations)[modify_pop_indices]*
                                                   dd_multipliers[modify_pop_indices])
        } else if (all(dispersal_target_n_k$threshold > dispersal_target_n_k$cutoff)) { # seek company /
          dd_multipliers[which(abundance_on_capacity <= dispersal_target_n_k$cutoff)] <- 0
          modify_pop_indices <- which(abundance_on_capacity < dispersal_target_n_k$threshold & dd_multipliers > 0)
          dd_multipliers[modify_pop_indices] <- ((abundance_on_capacity[modify_pop_indices] -
                                                    array(dispersal_target_n_k$cutoff, populations)[modify_pop_indices])/
                                                   array(dispersal_target_n_k$threshold - dispersal_target_n_k$cutoff,
                                                         populations)[modify_pop_indices]*
                                                   dd_multipliers[modify_pop_indices])
        }
      }

      # Select multipliers via target populations for non-zero occupied dispersals
      selected_dd_multipliers_list <- map(1:length(dispersal_target_pop_map_list), function(i) {
        dd_multipliers[dispersal_target_pop_map_list[[i]][, occupied_indices_list[[i]]][occupied_dispersal_indices_list[[i]]]]
      })

      # Apply modifying multipliers to dispersals
      modify_indices_list <- map(selected_dd_multipliers_list, \(x) which(x < 1))
      if (sum(lengths(modify_indices_list))) {
        modify_dispersal_indices_list <- map2(occupied_dispersal_indices_list, modify_indices_list,
                                              \(x, y) x[y])
        occupied_dispersals_list <- mapply(function(x, y, z, q) {
          x[y] <- x[y]*z[q]; x
        }, occupied_dispersals_list, modify_dispersal_indices_list,
        selected_dd_multipliers_list, modify_indices_list,
        SIMPLIFY = FALSE)
        occupied_dispersal_indices_list <- lapply(occupied_dispersals_list,
          function(x) which(as.logical(x))
        ) # > 0
      }

    } # dispersal depends on target pop N, K or N/K?

    for (segment in 1:(stages*compartments)) {

      if (!dispersal_stages_expanded[segment]) {
        next
      }

      # Disperser generation via abundance and corresponding dispersal rates
      occupied_abundance <- segment_abundance[segment, occupied_indices_list[[segment]]]
      occupied_abundance_rep <- segment_abundance[rep(segment, dispersal_compact_rows_list[[segment]]), occupied_indices_list[[segment]]]
      dispersers <- array(0, c(dispersal_compact_rows_list[[segment]], occupied_population_list[[segment]]))

      # Generate dispersers
      if (demographic_stochasticity) { # via binomial distribution
        dispersers[occupied_dispersal_indices_list[[segment]]] <- stats::rbinom(length(occupied_dispersal_indices_list[[segment]]),
                                                                                occupied_abundance_rep[occupied_dispersal_indices_list[[segment]]],
                                                                                occupied_dispersals_list[[segment]])
      } else { # deterministic
        dispersers[occupied_dispersal_indices_list[[segment]]] <- round(occupied_abundance_rep[occupied_dispersal_indices_list[[segment]]]*
                                                                      occupied_dispersals_list[[segment]])
      }

      # Calculate emigrants
      emigrants <- array(.colSums(dispersers, m = dispersal_compact_rows_list[[segment]], n = occupied_population_list[[segment]]))

      # Check consistency of emigrants (not to exceed abundances)
      excessive_indices <- which(emigrants > occupied_abundance)
      if (length(excessive_indices) > 0) { # reduce emigrants to equal abundance via random sampling
        for (excessive_index in excessive_indices) {
          excessive_rows <- which(as.logical(dispersers[, excessive_index])) # > 0
          excessive_dispersers <- dispersers[excessive_rows, excessive_index]
          disperser_reduction <- emigrants[excessive_index] - occupied_abundance[excessive_index]
          for (remove_row_index in rep(excessive_rows,
                                      times = excessive_dispersers)[sample(sum(excessive_dispersers),
                                                                            size = disperser_reduction)]) {
            dispersers[remove_row_index, excessive_index] <- dispersers[remove_row_index, excessive_index] - 1
          }
        }
        emigrants[excessive_indices] <- occupied_abundance[excessive_indices]
      }

      # Update occupied segment abundance
      segment_abundance[segment, occupied_indices_list[[segment]]] <- segment_abundance[segment, occupied_indices_list[[segment]]] - emigrants

      # Calculate immigrants via dispersal immigrant map
      disperser_indices <- which(as.logical(dispersers)) # > 0
      immigrant_array <- array(0, c(dispersal_compact_rows_list[[segment]], populations))
      immigrant_array[dispersal_immigrant_map_list[[segment]][, occupied_indices_list[[segment]]][disperser_indices]] <- dispersers[disperser_indices]
      immigrants <- .colSums(immigrant_array, m = dispersal_compact_rows_list[[segment]], n = populations)

      # Update population abundances
      segment_abundance[segment,] <- segment_abundance[segment,] + immigrants
    } # end dispersal

    # Perform additional dispersal for overcrowded cells (only to cells with room)
    if ((dispersal_depends_on_target_pop_n && all(dispersal_target_n$threshold < dispersal_target_n$cutoff)) ||
        (dispersal_depends_on_target_pop_n_k && all(dispersal_target_n_k$threshold < dispersal_target_n_k$cutoff))) {

      # Flags for dependencies
      depends_on_target_pop_n <- (dispersal_depends_on_target_pop_n && all(dispersal_target_n$threshold < dispersal_target_n$cutoff))
      depends_on_target_pop_n_k <- (dispersal_depends_on_target_pop_n_k && all(dispersal_target_n_k$threshold < dispersal_target_n_k$cutoff))

      # Get all updated dispersal rates
      dispersals <- pmap(list(rate = dispersal_compact_matrix_tm_list,
                              update = occupied_dispersals_list,
                              occupied = occupied_indices_list,
                              l = dispersal_stages_expanded),
                         \(rate, update, occupied, l) {
                           if (l) {
                             rate[, occupied] <- update
                           }
                           return(rate)
                         })

      # Identify overcrowded cells
      density_abundance <- .colSums(segment_abundance, m = stages*compartments, n = populations)
      if (depends_on_target_pop_n) {
        excessive_indices <- which(density_abundance > dispersal_target_n$cutoff)
      }
      if (depends_on_target_pop_n_k) {
        excessive_indices <- unique(c(excessive_indices,
                                      which(density_abundance/carrying_capacity > dispersal_target_n_k$cutoff)))
      }
      # Disperse excess from each overcrowded cell (in random order)
      for (segment in 1:(stages*compartments)) {

        if (!dispersal_stages_expanded[segment]) {
          next
        }

        excessive_indices_segment <- excessive_indices[segment_abundance[segment, excessive_indices] > 0]

        for (excessive_index in excessive_indices_segment[sample(length(excessive_indices_segment))]) {

          dispersal_indices <- which(dispersals[[segment]][, excessive_index] > 0)
          target_indices <- dispersal_target_pop_map_list[[segment]][, excessive_index][dispersal_indices]
          if (depends_on_target_pop_n && depends_on_target_pop_n_k) {
            indices_with_room <- which((density_abundance < dispersal_target_n$cutoff &
                                          (density_abundance + 1)/carrying_capacity <= dispersal_target_n_k$cutoff)[target_indices])
          } else if (depends_on_target_pop_n) {
            indices_with_room <- which((density_abundance < dispersal_target_n$cutoff)[target_indices])
          } else if (depends_on_target_pop_n_k) {
            indices_with_room <- which(((density_abundance + 1)/carrying_capacity <= dispersal_target_n_k$cutoff)[target_indices])
          }
          dispersal_indices <- dispersal_indices[indices_with_room]
          target_indices <- target_indices[indices_with_room]

          # Disperse excess one at a time sampled via the cell segment abundance distribution
          abundance_excess <- 0
          if (depends_on_target_pop_n) {
            abundance_excess <- density_abundance[excessive_index] - dispersal_target_n$cutoff
          }
          if (depends_on_target_pop_n_k) {
            abundance_excess <- max(abundance_excess, density_abundance[excessive_index] - floor(dispersal_target_n_k$cutoff*carrying_capacity[excessive_index]))
          }

          if (length(target_indices)) {

            # Sample target cell
            target_i <- target_indices[sample(length(target_indices), size = 1,
                                              prob = dispersals[[segment]][dispersal_indices, excessive_index])]

            # Perform dispersal
            segment_abundance[segment, excessive_index] <- segment_abundance[segment, excessive_index] - 1 # emigrant
            segment_abundance[segment, target_i] <- segment_abundance[segment, target_i] + 1 # immigrant

            # Update target density abundance and potential targets if it becomes full
            density_abundance[target_i] <- density_abundance[target_i] + 1
            if (((depends_on_target_pop_n &&
                  length(dispersal_target_n$cutoff) == 1 &&
                  density_abundance[target_i] >= dispersal_target_n$cutoff) ||
                 (depends_on_target_pop_n &&
                  length(dispersal_target_n$cutoff) > 1 &&
                  density_abundance[target_i] >= dispersal_target_n$cutoff[target_i])) ||
                ((depends_on_target_pop_n_k &&
                  density_abundance[target_i]/carrying_capacity[target_i] >= dispersal_target_n_k$cutoff))) {
              # remove from potential targets
              full_index <- which(target_indices == target_i)
              target_indices <- target_indices[-full_index]
              dispersal_indices <- dispersal_indices[-full_index]
            }
          }
        }
      }

    }

    return(segment_abundance)
  }

  return(dispersal_function)
}
