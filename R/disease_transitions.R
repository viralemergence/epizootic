#' Nested functions for stage- and compartment-based population transitions in
#' an outbreak.
#'
#' Modular functions for the disease simulator to transition populations
#' between stages and within disease compartments.
#'
#' @param stages Number of life cycle stages.
#' @param compartments Number of disease compartments.
#' @return Transition calculation function that takes as input:
#'  \describe{
#'     \item{\code{segment_abundance}}{Matrix of (current) abundance for each
#'     stage-compartment combo (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied
#'     at (current) time step.}
#'  }
#' @export disease_transitions

disease_transitions <- function(stages, compartments) {

  adult_indices <- seq(stages, stages*compartments, stages)
  non_adult_indices <- 1:(stages*compartments) |> _[-adult_indices]

  calculate <- function(segment_abundance, occupied_indices) {
    if (nrow(segment_abundance)==1) {
      new_matrix <- segment_abundance
    } else {
      new_matrix <- matrix(0, nrow = nrow(segment_abundance),
                           ncol = ncol(segment_abundance))
      new_matrix[non_adult_indices + 1, occupied_indices] <-
        segment_abundance[non_adult_indices, occupied_indices]
      new_matrix[adult_indices, occupied_indices] <-
        segment_abundance[adult_indices, occupied_indices] +
        new_matrix[adult_indices, occupied_indices]
    }
    return(new_matrix)
  }
  return(calculate)
}
