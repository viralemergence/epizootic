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

  # Convert the units of anything seasonal to now be daily
  fecundity <- replicate(active_pops, list(fecundity)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(fecundity_unit)] <- x[as.logical(fecundity_unit)] |> (`/`)(y)
      return(x)
    })
  mortality <- replicate(active_pops, list(mortality)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(mortality_unit)] <- x[as.logical(mortality_unit)] |> (`/`)(y)
      return(x)
    })
  transmission <- replicate(active_pops, list(transmission)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(transmission_unit)] <- x[as.logical(transmission_unit)] |>
        (`/`)(y)
      return(x)
    })
  recovery <- replicate(active_pops, list(recovery)) |>
    map2(breeding_season_length[occupied_indices], \(x, y) {
      x[as.logical(recovery_unit)] <- x[as.logical(recovery_unit)] |> (`/`)(y)
      return(x)
    })

  # Set up initial vectors
  population_list <- array_branch(segment_abundance[, occupied_indices], 2)

  population_new <- pmap(list(population_list, breeding_season_length,
                              mortality, transmission, recovery, fecundity,
                              abundance_threshold, carrying_capacity,
                              season = "breeding"), aspatial_siri)

  # Assign populations to occupied indices in segment_abundance
  for(i in occupied_indices) {
    segment_abundance[, i] <- population_new[[i]]
  }

  return(population_new)
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

  if (length(compare(names(initial_vector),
                     c("Sj", "Sa", "I1j", "I1a", "Rj", "Ra", "I2j", "I2a")))>0) {
    stop("Please label initial_vector correctly")
  }

  if (length(compare(names(parms),
                     c("season_length", "carrying_capacity", "beta_Sj", "beta_Sa",
                       "beta_Rj", "beta_Ra", "mortality_Sj", "mortality_Sa",
                       "mortality_I1j", "mortality_I1a", "mortality_I2j", "mortality_I2a",
                       "recovery_I1", "recovery_I2")))>0) {
    stop("Please label parms correctly")
  }

  if (anyNA(initial_vector)) {
    stop("initial_vector cannot have missing values")
  }

  if (anyNA(parms)) {
    stop("parms cannot have missing values")
  }

  # Unpack parameters
  winter_length <- parms["season_length"]
  beta_Sj <- parms["beta_Sj"]
  mSj <- parms["mortality_Sj"]
  beta_Sa <- parms["beta_Sa"]
  mSa <- parms["mortality_Sa"]
  mI1j <- parms["mortality_I1j"]
  recovery_I1 <- parms["recovery_I1"]
  mI1a <- parms["mortality_I1a"]
  recovery_I2 <- parms["recovery_I2"]
  beta_Rj <- parms["beta_Rj"]
  beta_Ra <- parms["beta_Ra"]
  mI2j <- parms["mortality_I2j"]
  mI2a <- parms["mortality_I2a"]
  K <- parms["carrying_capacity"]

  state_list <- c(list(initial_vector), replicate(winter_length, c(
    Sa = NA,
    Sj = NA,
    I1j = NA,
    I1a = NA,
    Rj = NA,
    Ra = NA,
    I2j = NA,
    I2a = NA
  ), simplify = F))

  for (t in 1:winter_length) {
    # Unpack states
    Sa <- state_list[[t]]["Sa"]
    Sj <- state_list[[t]]["Sj"]
    I1j <- state_list[[t]]["I1j"]
    I1a <- state_list[[t]]["I1a"]
    Rj <- state_list[[t]]["Rj"]
    Ra <- state_list[[t]]["Ra"]
    I2j <- state_list[[t]]["I2j"]
    I2a <- state_list[[t]]["I2a"]

    # define N
    N <- sum(Sj + Sa + I1j + I1a + Rj + Ra + I2j + I2a, na.rm = T)

    if (N > K) {
      N <- K
    }

    Sa_mortality <- (1+N/K)*mSa
    Sj_mortality <- (1+N/K)*mSj

    if (Sj_mortality > 1) {
      Sj_mortality <- 1
    }

    if (Sa_mortality > 1) {
      Sa_mortality <- 1
    }

    infection1_juv <- rbinom(1, prob = beta_Sj, Sj*(I1j + I2j + I1a + I2a))
    infection1_juv <- if_else(infection1_juv > Sj, Sj, infection1_juv)
    infection1_adult <- rbinom(1, prob = beta_Sa, Sa*(I1j + I2j + I1a + I2a))
    infection1_adult <- if_else(infection1_adult > Sa, Sa, infection1_adult)
    recovery1_juv <- rbinom(1, prob = recovery_I1, I1j)
    recovery1_adult <- rbinom(1, prob = recovery_I1, I1a)
    recovery2_juv <- rbinom(1, prob = recovery_I2, I2j)
    recovery2_adult <- rbinom(1, prob = recovery_I2, I2a)
    infection2_juv <- rbinom(1, prob = beta_Rj, Rj*(I1j + I2j + I1a + I2a))
    infection2_juv <- if_else(infection2_juv > Rj, Rj, infection2_juv)
    infection2_adult <- rbinom(1, prob = beta_Ra, Ra*(I1j + I2j + I1a + I2a))
    infection2_adult <- if_else(infection2_adult > Ra, Ra, infection2_adult)
    susceptible_adult_death <- rbinom(1, prob = Sa_mortality, Sa - infection1_adult)
    recovered_juvenile_death <- rbinom(1, prob = Sj_mortality, Rj + recovery1_juv + recovery2_juv - infection2_juv)
    recovered_adult_death <- rbinom(1, prob = Sa_mortality, Ra + recovery1_adult + recovery2_adult - infection2_adult)
    susceptible_juvenile_death <- rbinom(1, prob = Sj_mortality, Sj - infection1_juv)
    infected1_juvenile_death <- rbinom(1, prob = (1+N/K)*mI1j, I1j + infection1_juv - recovery1_juv)
    infected1_adult_death <- rbinom(1, prob = (1+N/K)*mI1a, I1a + infection1_adult  - recovery1_adult)
    infected2_juvenile_death <- rbinom(1, prob = (1+N/K)*mI2j, I2j + infection2_juv - recovery2_juv)
    infected2_adult_death <- rbinom(1, prob = (1+N/K)*mI2a, I2a + infection2_adult - recovery2_adult)

    # Simulate with demographic stochasticity
    state_list[[t+1]]["Sj"] <- Sj - infection1_juv - susceptible_juvenile_death
    state_list[[t+1]]["Sa"] <- Sa - infection1_adult - susceptible_adult_death
    state_list[[t+1]]["I1j"] <- I1j + infection1_juv - recovery1_juv - infected1_juvenile_death
    state_list[[t+1]]["I1a"] <- I1a + infection1_adult  - recovery1_adult - infected1_adult_death
    state_list[[t+1]]["Rj"] <- Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death
    state_list[[t+1]]["Ra"] <- Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death
    state_list[[t+1]]["I2j"] <- I2j + infection2_juv - recovery2_juv - infected2_juvenile_death
    state_list[[t+1]]["I2a"] <- I2a + infection2_adult - recovery2_adult - infected2_adult_death
  }

  df <- state_list %>% bind_rows() %>% rowid_to_column("Day")

  return(df)
}

#' Helper Function for Seasonal SIRI Simulation
#'
#' This function is an internal one that does the aspatial simulations within
#' one population for one timestep, for any given season.
#'
#' @param initial_pop A vector of length 8 showing the initial abundance for
#' each combination of stage and compartment.
#' @param season_length The length of the season in days.
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

aspatial_siri <- function(initial_pop, season_length, mortality, transmission,
                          recovery, fecundity = NULL, abundance_threshold,
                          carrying_capacity,
                          season = c("breeding", "non-breeding")) {

  state <- matrix(0, nrow = 8, ncol = season_length + 1)
  state[,1] <- initial_pop

  for (t in 1:season_length) {
    # Unpack states
    Sa <- state[2, t]
    Sj <- state[1, t]
    I1j <- state[3, t]
    I1a <- state[4, t]
    Rj <- state[5, t]
    Ra <- state[6, t]
    I2j <- state[7, t]
    I2a <- state[8, t]

    N <- pmin(sum(state[, t]), carrying_capacity)

    if (N < abundance_threshold) {
      # Fill the rest of the population's state vectors with 0 and skip to the
      # end of the loop
      state[, t:(ncol(state))] <- 0
      break
    }

    dd_mortality <- (1+N/carrying_capacity) |> (`*`)(mortality) |>
      map_dbl(\(x) pmin(x, 1))

    # Simulate with demographic stochasticity
    if (season == "breeding") {
      births <- rpois(sum(Sa + I1a + Ra + I2a),
                      fecundity * (1 - N / carrying_capacity))
    }
    infection1_juv <- rbinom(1, prob = transmission[1],
                             Sj*(I1j + I2j + I1a + I2a))
    infection1_juv <- pmin(infection1_juv, Sj)
    infection1_adult <- rbinom(1, prob = transmission[2],
                               Sa*(I1j + I2j + I1a + I2a))
    infection1_adult <- pmin(infection1_adult, Sa)
    recovery1_juv <- rbinom(1, prob = recovery[1],
                            I1j)
    recovery1_adult <- rbinom(1, prob = recovery[2],
                              I1a)
    recovery2_juv <- rbinom(1, prob = recovery[3],
                            I2j)
    recovery2_adult <- rbinom(1, prob = recovery[4],
                              I2a)
    infection2_juv <- rbinom(1, prob = transmission[3],
                             Rj*(I1j + I2j + I1a + I2a))
    infection2_juv <- pmin(infection2_juv, Rj)
    infection2_adult <- rbinom(1, prob = transmission[4],
                               Ra*(I1j + I2j + I1a + I2a))
    infection2_adult <- pmin(infection2_adult, Ra)
    susceptible_adult_death <- rbinom(1, prob = dd_mortality[2],
                                      Sa - infection1_adult)
    recovered_juvenile_death <- rbinom(1, prob = dd_mortality[5],
                                       Rj + recovery1_juv + recovery2_juv -
                                         infection2_juv)
    recovered_adult_death <- rbinom(1, prob = dd_mortality[6],
                                    Ra + recovery1_adult + recovery2_adult -
                                      infection2_adult)
    susceptible_juvenile_death <- rbinom(1, prob = dd_mortality[1],
                                         Sj + sum(births) - infection1_juv)
    infected1_juvenile_death <- rbinom(1, prob = dd_mortality[3],
                                       I1j + infection1_juv - recovery1_juv)
    infected1_adult_death <- rbinom(1, prob = dd_mortality[4],
                                    I1a + infection1_adult  - recovery1_adult)
    infected2_juvenile_death <- rbinom(1, prob = dd_mortality[7],
                                       I2j + infection2_juv - recovery2_juv)
    infected2_adult_death <- rbinom(1, prob = dd_mortality[8],
                                    I2a + infection2_adult - recovery2_adult)

    # Record new states
    state[1, t + 1] <- if_else(season == "breeding",
                               Sj + sum(births) - infection1_juv -
                                 susceptible_juvenile_death,
                               Sj - infection1_juv - susceptible_juvenile_death)
    state[2, t + 1] <- Sa - infection1_adult - susceptible_adult_death
    state[3, t + 1] <- I1j + infection1_juv - recovery1_juv -
      infected1_juvenile_death
    state[4, t + 1] <- I1a + infection1_adult  - recovery1_adult -
      infected1_adult_death
    state[5, t + 1] <- Rj + recovery1_juv + recovery2_juv - infection2_juv -
      recovered_juvenile_death
    state[6, t + 1] <- Ra + recovery1_adult + recovery2_adult -
      infection2_adult - recovered_adult_death
    state[7, t + 1] <- I2j + infection2_juv - recovery2_juv -
      infected2_juvenile_death
    state[8, t + 1] <- I2a + infection2_adult - recovery2_adult -
      infected2_adult_death
  }
  return(state[, season_length + 1])
}
