#' Simulate a *Mycoplasma gallisepticum* Outbreak in the Breeding Season
#'
#' Simulate a *Mycoplasma gallisepticum* outbreak during the breeding season
#' day-by-day in a population of house finches (*Haemorhous mexicanus*). Uses a
#' SIRI model (Susceptible-Infected 1- Recovered-Infected 2+) and includes
#' demographic stochasticity in fecundity, mortality, and infection.
#'
#' The function can also handle the case in which there are no infected
#' individuals. The principal difference between this function and the one for
#' simulating an outbreak in the non-breeding season is that this one includes
#' fecundity.
#'
#' @importFrom stats rbinom rpois
#' @import dplyr
#' @import tibble
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
#' @return An abundance matrix with `populations` columns and
#' `stages*compartments` rows, updated from the `segment_abundance` input in
#' `inputs` according to demography and disease dynamics.
#' @export

siri_model_summer <- function(inputs) {
  inputs <- check_aspatial_siri_inputs(inputs)
  list2env(inputs, environment())
  active_pops <- length(occupied_indices)

  # Convert the units of anything seasonal to now be daily
  fecundity <- replicate(active_pops, list(fecundity)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(fecundity_unit)] <- x[as.logical(fecundity_unit)] |>
        (`/`)(y)
      return(x)
    })
  mortality <- replicate(active_pops, list(mortality)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(mortality_unit)] <- x[as.logical(mortality_unit)] |>
        (`/`)(y) |> pmin(1)
      return(x)
    })
  transmission <- replicate(active_pops, list(transmission)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(transmission_unit)] <- x[as.logical(transmission_unit)] |>
        (`/`)(y) |> pmin(1)
      return(x)
    })
  recovery <- replicate(active_pops, list(recovery)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(recovery_unit)] <- x[as.logical(recovery_unit)] |>
        (`/`)(y) |> pmin(1)
      return(x)
    })

  # Set up initial vectors
  if (active_pops > 1) {
    population_list <- array_branch(segment_abundance[, occupied_indices], 2)
  }
  if (active_pops == 1) {
    population_list <- list(segment_abundance[, occupied_indices])
  }
  if (active_pops == 0) {
    return(segment_abundance)
  }
  carrying_capacity <- carrying_capacity[occupied_indices]
  season_length <- breeding_season_length[occupied_indices]

  population_new <- pmap(list(initial_pop = population_list,
                              mortality = mortality,
                              transmission = transmission,
                              recovery = recovery,
                              fecundity = fecundity,
                              abundance_threshold = abundance_threshold,
                              carrying_capacity = carrying_capacity,
                              season_length = season_length,
                              season = "breeding"), aspatial_siri)

  # Assign populations to occupied indices in segment_abundance
  for(i in 1:length(occupied_indices)) {
    segment_abundance[, occupied_indices[i]] <- population_new[[i]]
  }

  return(segment_abundance)
}

#' Simulate a *Mycoplasma gallisepticum* Outbreak in the Non-Breeding Season
#'
#' Simulate a *Mycoplasma gallisepticum* outbreak during the non-breeding season day-by-day in a
#' population of house finches (*Haemorhous mexicanus*). Uses a SIRI model (Susceptible-Infected 1-
#' Recovered-Infected 2+) and includes demographic stochasticity in fecundity, mortality, and
#' infection.
#'
#' The function can also handle the case in which there are no infected individuals. The principal
#' difference between this function and the one for simulating an outbreak in the breeding season
#' is that this one does not include fecundity.
#'
#' @param inputs A nested list with named elements:
#'   \describe{
#'     \item{\code{replicates}}{Number of replicate simulation runs.}
#'     \item{\code{time_steps}}{Number of simulation time steps.}
#'     \item{\code{populations}}{Number of populations.}
#'     \item{\code{stages}}{Number of life cycle stages.}
#'     \item{\code{compartments}}{Number of disease compartments (e.g., 3 for a
#'     SIR model).}
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
#'     \item{\code{simulator}}{[`poems::SimulatorReference`] object with
#'     dynamically accessible \emph{attached} and \emph{results} lists.}
#'     \item{\code{additional attributes}}{Additional attributes when the
#'     transformation is optionally nested in a list.}
#'   }
#' @return An abundance matrix with `populations` columns and
#' `stages*compartments` rows, updated from the `segment_abundance` input in
#' `inputs` according to demography and disease dynamics.
#' @export

siri_model_winter <- function(inputs) {
  inputs <- check_aspatial_siri_inputs(inputs)
  list2env(inputs, environment())
  active_pops <- length(occupied_indices)
  season_length <- (365 - breeding_season_length)[occupied_indices]

  # Convert the units of anything seasonal to now be daily
  fecundity <- replicate(active_pops, list(fecundity)) |>
    map2(season_length, \(x, y) {
      x[as.logical(fecundity_unit)] <- x[as.logical(fecundity_unit)] |>
        (`/`)(y)
      return(x)
    })
  mortality <- replicate(active_pops, list(mortality)) |>
    map2(season_length, \(x, y) {
      x[as.logical(mortality_unit)] <- x[as.logical(mortality_unit)] |>
        (`/`)(y) |> pmin(1)
      return(x)
    })
  transmission <- replicate(active_pops, list(transmission)) |>
    map2(season_length, \(x, y) {
      x[as.logical(transmission_unit)] <- x[as.logical(transmission_unit)] |>
        (`/`)(y) |> pmin(1)
      return(x)
    })
  recovery <- replicate(active_pops, list(recovery)) |>
    map2(season_length, \(x, y) {
      x[as.logical(recovery_unit)] <- x[as.logical(recovery_unit)] |>
        (`/`)(y) |> pmin(1)
      return(x)
    })

  # Set up initial vectors
  if (active_pops > 1) {
    population_list <- array_branch(segment_abundance[, occupied_indices], 2)
  }
  if (active_pops == 1) {
    population_list <- list(segment_abundance[, occupied_indices])
  }
  if (active_pops == 0) {
    return(segment_abundance)
  }
  carrying_capacity <- carrying_capacity[occupied_indices]

  population_new <- pmap(list(initial_pop = population_list,
                              mortality = mortality,
                              transmission = transmission,
                              recovery = recovery,
                              fecundity = fecundity,
                              abundance_threshold = abundance_threshold,
                              carrying_capacity = carrying_capacity,
                              season_length = season_length,
                              season = "non-breeding"), aspatial_siri)

  # Assign populations to occupied indices in segment_abundance
  for(i in 1:length(occupied_indices)) {
    segment_abundance[, occupied_indices[i]] <- population_new[[i]]
  }

  return(segment_abundance)
}
