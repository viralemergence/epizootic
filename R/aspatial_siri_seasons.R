#' Simulate a *Mycoplasma gallisepticum* Outbreak in the Breeding Season
#'
#' Simulate a *Mycoplasma gallisepticum* outbreak during the breeding season day-by-day in a
#' population of house finches (*Haemorhous mexicanus*). Uses a SIRI model (Susceptible-Infected 1-
#' Recovered-Infected 2+) and includes demographic stochasticity in fecundity, mortality, and
#' infection.
#'
#' The function can also handle the case in which there are no infected individuals. The principal
#' difference between this function and the one for simulating an outbreak in the non-breeding season
#' is that this one includes fecundity.
#'
#' @import waldo
#' @import dplyr
#' @param initial_vector Represents the state of the population at the start of the breeding season.
#' A list with named elements, each one a numeric of length 1:
#' \describe{
#'     \item{\code{Sj}}{Susceptible juveniles.}
#'     \item{\code{Sa}}{Susceptible adults.}
#'     \item{\code{I1j}}{Juveniles infected for the first time.}
#'     \item{\code{I1a}}{Adults infected for the first time.}
#'     \item{\code{Rj}}{Recovered juveniles.}
#'     \item{\code{Ra}}{Recovered adults.}
#'     \item{\code{I2j}}{Juveniles infected for the second+ time.}
#'     \item{\code{I2a}}{Adults infected for the second+ time.}
#' }
#' @param parms A vector of simulation parameters:
#' \describe{
#'     \item{\code{season_length}}{Length of the breeding season in days.}
#'     \item{\code{carrying_capacity}}{Carrying capacity of finches in the population.}
#'     \item{\code{birth}}{Daily fecundity.}
#'     \item{\code{beta_Sj}}{Daily transmission for susceptible juveniles in breeding season.}
#'     \item{\code{beta_Sa}}{Daily transmission for susceptible adults in breeding season.}
#'     \item{\code{beta_Rj}}{Daily transmission for recovered juveniles in breeding season.}
#'     \item{\code{beta_Ra}}{Daily transmission for recovered adults in breeding season.}
#'     \item{\code{mortality_Sj}}{Daily mortality of susceptible juveniles in breeding season.}
#'     \item{\code{mortality_Sa}}{Daily mortality of susceptible adults in breeding season.}
#'     \item{\code{mortality_I1j}}{Daily mortality of juveniles infected for the first time.}
#'     \item{\code{mortality_I1a}}{Daily mortality of adults infected for the first time.}
#'     \item{\code{mortality_I2j}}{Daily mortality of juveniles infected for the second+ time.}
#'     \item{\code{mortality_I2a}}{Daily mortality of adults infected for the second+ time.}
#'     \item{\code{recovery_I1}}{Daily recovery rate of birds infected for the first time.}
#'     \item{\code{recovery_I2}}{Daily recover rate of birds infected for the second+ time.}
#' }
#' @param ... Does nothing. A placeholder for future code improvements.
#' @return A data frame where each row is a day of the simulation and each column is one of the
#' population states detailed in `initial_vector`.

siri_model_summer <- function(initial_vector, parms, ...) {
  if (length(compare(names(initial_vector),
                     c("Sj", "Sa", "I1j", "I1a", "Rj", "Ra", "I2j", "I2a")))>0) {
    stop("Please label initial_vector correctly")
  }

  if (length(compare(names(parms),
                     c("season_length", "carrying_capacity", "birth", "beta_Sj", "beta_Sa",
                       "beta_Rj", "beta_Ra", "mortality_Sj", "mortality_Sa",
                       "mortality_I1j", "mortality_I1a", "mortality_I2j", "mortality_I2a",
                       "recovery_I1", "recovery_I2")))>0) {
    stop("Please label parms correctly")
  }

  if (any(is.na(initial_vector))) {
    stop("initial_vector cannot have missing values")
  }

  if (any(is.na(parms))) {
    stop("parms cannot have missing values")
  }

  # Unpack parameters
  summer_length <- parms["season_length"]
  daily_egg <- parms["birth"]
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

  state_list <- c(list(y), replicate(i, c(
    Sa = NA,
    Sj = NA,
    I1j = NA,
    I1a = NA,
    Rj = NA,
    Ra = NA,
    I2j = NA,
    I2a = NA
  ), simplify = F))

  for (t in 1:i) {
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

    Sj_mortality <- (1+N/K)*mSj
    Sa_mortality <- (1+N/K)*mSa

    if (Sj_mortality > 1) {
      Sj_mortality <- 1
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
    births <- try(rpois(sum(Sa + I1a + Ra + I2a), daily_egg*(1-N/K)), silent = TRUE)
    susceptible_juvenile_death <- rbinom(1, prob = Sj_mortality, Sj + sum(births) - infection1_juv)
    infected1_juvenile_death <- rbinom(1, prob = (1+N/K)*mI1j, I1j + infection1_juv - recovery1_juv)
    infected1_adult_death <- rbinom(1, prob = (1+N/K)*mI1a, I1a + infection1_adult  - recovery1_adult)
    infected2_juvenile_death <- rbinom(1, prob = (1+N/K)*mI2j, I2j + infection2_juv - recovery2_juv)
    infected2_adult_death <- rbinom(1, prob = (1+N/K)*mI2a, I2a + infection2_adult - recovery2_adult)

    if ("try-error" %in% class(births)) {
      break
    } else {

      # Simulate with demographic stochasticity
      state_list[[t+1]]["Sj"] <- Sj + sum(births) - infection1_juv - susceptible_juvenile_death
      state_list[[t+1]]["Sa"] <- Sa - infection1_adult - susceptible_adult_death
      state_list[[t+1]]["I1j"] <- I1j + infection1_juv - recovery1_juv - infected1_juvenile_death
      state_list[[t+1]]["I1a"] <- I1a + infection1_adult  - recovery1_adult - infected1_adult_death
      state_list[[t+1]]["Rj"] <- Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death
      state_list[[t+1]]["Ra"] <- Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death
      state_list[[t+1]]["I2j"] <- I2j + infection2_juv - recovery2_juv - infected2_juvenile_death
      state_list[[t+1]]["I2a"] <- I2a + infection2_adult - recovery2_adult - infected2_adult_death
    }
  }

  df <- state_list %>% bind_rows() %>% rowid_to_column("Day")

  return(df)
}
