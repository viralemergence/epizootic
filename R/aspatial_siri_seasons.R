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
#'     \item{\code{simulator}}{\code{\link{SimulatorReference}} object with
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

  # Convert the units of anything daily to now be seasonal
  fecundity <- replicate(active_pops, list(fecundity)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[!as.logical(fecundity_unit)] <- x[!as.logical(fecundity_unit)] |>
        (`*`)(y)
      return(x)
    })
  mortality <- replicate(active_pops, list(mortality)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[!as.logical(mortality_unit)] <- x[!as.logical(mortality_unit)] |>
        (`*`)(y) |> pmin(1)
      return(x)
    })
  transmission <- replicate(active_pops, list(transmission)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[!as.logical(transmission_unit)] <- x[!as.logical(transmission_unit)] |>
        (`*`)(y) |> pmin(1)
      return(x)
    })
  recovery <- replicate(active_pops, list(recovery)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[!as.logical(recovery_unit)] <- x[!as.logical(recovery_unit)] |>
        (`*`)(y) |> pmin(1)
      return(x)
    })

  # Set up initial vectors
  population_list <- array_branch(segment_abundance[, occupied_indices], 2)
  carrying_capacity <- carrying_capacity[occupied_indices]

  population_new <- pmap(list(initial_pop = population_list,
                              mortality = mortality,
                              transmission = transmission,
                              recovery = recovery,
                              fecundity = fecundity,
                              abundance_threshold = abundance_threshold,
                              carrying_capacity = carrying_capacity,
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
#'     \item{\code{simulator}}{\code{\link{SimulatorReference}} object with
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
  season_length <- 365 - breeding_season_length

  # Convert the units of anything daily to now be seasonal
  fecundity <- replicate(active_pops, list(fecundity)) |>
    map2(season_length[occupied_indices], \(x, y) {
      x[!as.logical(fecundity_unit)] <- x[!as.logical(fecundity_unit)] |>
        (`*`)(y)
      return(x)
    })
  mortality <- replicate(active_pops, list(mortality)) |>
    map2(season_length[occupied_indices], \(x, y) {
      x[!as.logical(mortality_unit)] <- x[!as.logical(mortality_unit)] |>
        (`*`)(y) |> pmin(1)
      return(x)
    })
  transmission <- replicate(active_pops, list(transmission)) |>
    map2(season_length[occupied_indices], \(x, y) {
      x[!as.logical(transmission_unit)] <- x[!as.logical(transmission_unit)] |>
        (`*`)(y) |> pmin(1)
      return(x)
    })
  recovery <- replicate(active_pops, list(recovery)) |>
    map2(season_length[occupied_indices], \(x, y) {
      x[!as.logical(recovery_unit)] <- x[!as.logical(recovery_unit)] |>
        (`*`)(y) |> pmin(1)
      return(x)
    })

  # Set up initial vectors
  population_list <- array_branch(segment_abundance[, occupied_indices], 2)
  carrying_capacity <- carrying_capacity[occupied_indices]

  population_new <- pmap(list(initial_pop = population_list,
                              mortality = mortality,
                              transmission = transmission,
                              recovery = recovery,
                              fecundity = fecundity,
                              abundance_threshold = abundance_threshold,
                              carrying_capacity = carrying_capacity,
                              season = "non-breeding"), aspatial_siri)

  # Assign populations to occupied indices in segment_abundance
  for(i in 1:length(occupied_indices)) {
    segment_abundance[, occupied_indices[i]] <- population_new[[i]]
  }

  return(segment_abundance)
}

#' Helper function: SIRI seasonal simulator
#'
#' @param initial_pop A vector of length 8 showing the initial abundance for
#' each combination of stage and compartment.
#' @param mortality A vector of length 8 with the mortality rates for each
#' stage and compartment in the season in question.
#' @param transmission A vector of length 4 with the transmission rates for each
#' susceptible/recovered stage in the season in question.
#' @param recovery A vector of length 4 with the transmission rates for each
#' susceptible/recovered stage in the season in question.
#' @param fecundity Only necessary when `season = "breeding"` (see below).
#' Default NULL. A single numeric with the daily fecundity of adults.
#' @param abundance_threshold A quasi-extinction threshold below which a
#' population becomes extinct.
#' @param carrying_capacity A single numeric that indicates the carrying
#' capacity of the population in this season.
#' @param season Either "breeding" or "non-breeding."
#' @return A vector of length 8 showing the abundance for each combination of
#' stage and compartment at the end of the season.
#' @export

aspatial_siri <- function(initial_pop,
                          carrying_capacity,
                          abundance_threshold,
                          mortality,
                          transmission,
                          recovery,
                          fecundity,
                          season) {
  # Unpack initial_pops
  Sa <- initial_pop[2]
  Sj <- initial_pop[1]
  I1j <- initial_pop[3]
  I1a <- initial_pop[4]
  Rj <- initial_pop[5]
  Ra <- initial_pop[6]
  I2j <- initial_pop[7]
  I2a <- initial_pop[8]

  total_pop <- sum(initial_pop)

  N <- pmin(total_pop, carrying_capacity)

  if (N < abundance_threshold) {
    # Fill everything with 0
    final_state <- rep(0, 8)
  } else {

  dd_multiplier <- 1.0 + N / carrying_capacity

  dd_mortality <- mortality * pmin(dd_multiplier, 1.0)

  new_juv <- 0.0

  if (season == "breeding") {
    adults <- Sa + I1a + Ra + I2a
    births <- rpois(adults, fecundity * (1.0 - N / carrying_capacity))
    new_juv <- sum(births)
  }

  infection1_juv_raw <- rbinom(1, Sj * (I1j + I2j + I1a + I2a), transmission[1])
  infection1_juv <- pmin(infection1_juv_raw, Sj)

  infection1_adult_raw <- rbinom(1, Sa * (I1j + I2j + I1a + I2a), transmission[2])
  infection1_adult <- pmin(infection1_adult_raw, Sa)

  susceptible_adult_death <- rbinom(1, Sa - infection1_adult, dd_mortality[2])

  if (season == "breeding") {
    susceptible_juvenile_death <- rbinom(1, Sj + new_juv - infection1_juv, dd_mortality[1])
  } else {
    susceptible_juvenile_death <- rbinom(1, Sj - infection1_juv, dd_mortality[1])
  }

  infected1_juvenile_death <- rbinom(1, I1j + infection1_juv, dd_mortality[3])

  infected1_adult_death <- rbinom(1, I1a + infection1_adult, dd_mortality[4])

  recovery1_juv_raw <- rbinom(1, I1j + infection1_juv - infected1_juvenile_death, recovery[3])
  recovery1_juv <- pmin(recovery1_juv_raw, I1j + infection1_juv - infected1_juvenile_death)

  recovery1_adult_raw <- rbinom(1, I1a + infection1_adult - infected1_adult_death, recovery[4])
  recovery1_adult <- pmin(recovery1_adult_raw, I1a + infection1_adult - infected1_adult_death)

  infection2_juv_raw <- rbinom(1, (Rj + recovery1_juv) * (I1j + I2j + I1a + I2a), transmission[3])
  infection2_juv <- pmin(infection2_juv_raw, Rj + infection1_juv)

  infection2_adult_raw <- rbinom(1, (Ra + recovery1_adult) * (I1j + I2j + I1a + I2a), transmission[4])
  infection2_adult <- pmin(infection2_adult_raw, Ra + recovery1_adult)

  infected2_juvenile_death <- rbinom(1, I2j + infection2_juv, dd_mortality[7])

  infected2_adult_death <- rbinom(1, I2a + infection2_adult, dd_mortality[8])

  recovery2_juv_raw <- rbinom(1, I2j + infection2_juv - infected2_juvenile_death, recovery[3])
  recovery2_juv <- pmin(recovery2_juv_raw, I2j + infection2_juv - infected2_juvenile_death)

  recovery2_adult_raw <- rbinom(1, I2a + infection2_adult - infected2_adult_death, recovery[4])
  recovery2_adult <- pmin(recovery2_adult_raw, I2a + infection2_adult - infected2_adult_death)

  recovered_juvenile_death <- rbinom(1, Rj + recovery1_juv + recovery2_juv - infection2_juv, dd_mortality[5])

  recovered_adult_death <- rbinom(1, Ra + recovery1_adult + recovery2_adult - infection2_adult, dd_mortality[6])

  # Update state for the next time step
  if (season == "breeding") {
    final_state <- c(Sj + new_juv - infection1_juv - susceptible_juvenile_death,
                     Sa - infection1_adult - susceptible_adult_death,
                     I1j + infection1_juv - recovery1_juv - infected1_juvenile_death,
                     I1a + infection1_adult - recovery1_adult - infected1_adult_death,
                     Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death,
                     Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death,
                     I2j + infection2_juv - recovery2_juv - infected2_juvenile_death,
                     I2a + infection2_adult - recovery2_adult - infected2_adult_death)
  } else {
    final_state <- c(Sj - infection1_juv - susceptible_juvenile_death,
                     Sa - infection1_adult - susceptible_adult_death,
                     I1j + infection1_juv - recovery1_juv - infected1_juvenile_death,
                     I1a + infection1_adult - recovery1_adult - infected1_adult_death,
                     Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death,
                     Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death,
                     I2j + infection2_juv - recovery2_juv - infected2_juvenile_death,
                     I2a + infection2_adult - recovery2_adult - infected2_adult_death)
  }
  }

  return(final_state)
}
