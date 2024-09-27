#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Helper Function for Seasonal SIRI Simulation
//'
//' This function is an internal one that does the aspatial simulations within
//' one population for one timestep, for any given season.
//'
//' @name aspatial_siri
//'
//' @param initial_pop A vector of length 8 showing the initial abundance for
//' each combination of stage and compartment.
//' @param season_length The length of the season in days.
//' @param mortality A vector of length 8 with the mortality rates for each
//' stage and compartment in the season in question.
//' @param transmission A vector of length 8 with the transmission rates for each
//' stage in the season in question.
//' @param recovery A vector of length 8 with the recovery rates for each
//' infected stage in the season in question.
//' @param fecundity A vector of length 8 with the fecundity for each
//' reproductive segment.
//' @param abundance_threshold A quasi-extinction threshold below which a
//' population becomes extinct.
//' @param carrying_capacity A single numeric that indicates the carrying
//' capacity of the population in this season.
//' @param season Either "breeding" or "non-breeding."
//' @return A vector of length 8 showing the abundance for each combination of
//' stage and compartment at the end of the season.
//'
//' @examples
//' aspatial_siri(
//'  initial_pop = c(50000, 50000, 0, 1, 0, 0, 0, 0),
//'  season_length = 100,
//'  mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04),
//'  fecundity = c(0, 15/182, 0, 15/182, 0, 15/182, 0, 15/182),
//'  transmission = c(0.00002, 0.00001, 0, 0, 7.84e-06, 3.92e-06, 0, 0),
//'  recovery = c(0, 0, 0.05714286, 0.05714286, 0, 0, 0.1, 0.1),
//'  carrying_capacity = 150000,
//'  abundance_threshold = 10,
//'  season = "breeding"
//' )
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector aspatial_siri(Rcpp::NumericVector initial_pop,
                                  int season_length,
                                  Rcpp::NumericVector mortality,
                                  Rcpp::NumericVector transmission,
                                  Rcpp::NumericVector recovery,
                                  Rcpp::NumericVector fecundity,
                                  double abundance_threshold,
                                  double carrying_capacity,
                                  const std::string& season) {

  int n_stages = 8;
  arma::mat state(n_stages, season_length + 1);
  Rcpp::NumericVector dd_mortality(n_stages);
  bool isBreedingSeason = (season == "breeding");
  double birth_rate = fecundity[1];

  state.col(0) = Rcpp::as<arma::vec>(initial_pop);
  for (int t = 0; t < season_length; t++) {
    // Unpack states
    double Sj = state(0, t);
    double Sa = state(1, t);
    double I1j = state(2, t);
    double I1a = state(3, t);
    double Rj = state(4, t);
    double Ra = state(5, t);
    double I2j = state(6, t);
    double I2a = state(7, t);

    double N = std::min(arma::accu(state.col(t)), carrying_capacity);

    if (N < abundance_threshold) {
      state.cols(t + 1, season_length).zeros();
      break;
    }

    for (int i = 0; i < n_stages; i++) {
      dd_mortality[i] = std::min((1.0 + N / carrying_capacity) * mortality[i], 1.0);
    }

    double new_juv = 0;

    if (isBreedingSeason) {
      double adults = Sa + I1a + Ra + I2a;
      double birth_rate_adj = birth_rate * (1.0 - N / carrying_capacity);
      new_juv = Rcpp::rpois(1, adults * birth_rate_adj)[0];
    }

    bool hasInfections = (I1j + I1a + I2j + I2a) > 0;

    if (!hasInfections) {
      // Only update Sj and Sa, skipping all infection and recovery calculations
      state(0, t + 1) = isBreedingSeason ? Sj + new_juv - Rcpp::rbinom(1, Sj + new_juv, dd_mortality[0])[0] : Sj - Rcpp::rbinom(1, Sj, dd_mortality[0])[0];
      state(1, t + 1) = Sa - Rcpp::rbinom(1, Sa, dd_mortality[1])[0];
      state(2, t + 1) = I1j;
      state(3, t + 1) = I1a;
      state(4, t + 1) = Rj - Rcpp::rbinom(1, Rj, dd_mortality[4])[0];
      state(5, t + 1) = Ra - Rcpp::rbinom(1, Ra, dd_mortality[5])[0];
      state(6, t + 1) = I2j;
      state(7, t + 1) = I2a;
      continue;
    }

    double infection1_juv = std::min(Rcpp::rbinom(1, Sj * (I1j + I2j + I1a + I2a), transmission[0])[0], Sj);
    double infection1_adult = std::min(Rcpp::rbinom(1, Sa * (I1j + I2j + I1a + I2a), transmission[1])[0], Sa);
    double susceptible_adult_death = Rcpp::rbinom(1, Sa - infection1_adult, dd_mortality[1])[0];
    double susceptible_juvenile_death = Rcpp::rbinom(1, isBreedingSeason ? Sj + new_juv - infection1_juv : Sj - infection1_juv, dd_mortality[0])[0];
    double infected1_juvenile_death = Rcpp::rbinom(1, I1j + infection1_juv, dd_mortality[2])[0];
    double infected1_adult_death = Rcpp::rbinom(1, I1a + infection1_adult, dd_mortality[3])[0];
    double recovery1_juv = std::min(Rcpp::rbinom(1, I1j + infection1_juv - infected1_juvenile_death, recovery[2])[0], I1j + infection1_juv - infected1_juvenile_death);
    double recovery1_adult = std::min(Rcpp::rbinom(1, I1a + infection1_adult - infected1_adult_death, recovery[3])[0], I1a + infection1_adult - infected1_adult_death);
    double infection2_juv = std::min(Rcpp::rbinom(1, (Rj + recovery1_juv) * (I1j + I2j + I1a + I2a), transmission[4])[0], Rj + recovery1_juv);
    double infection2_adult = std::min(Rcpp::rbinom(1, (Ra + recovery1_adult) * (I1j + I2j + I1a + I2a), transmission[5])[0], Ra + recovery1_adult);
    double infected2_juvenile_death = Rcpp::rbinom(1, I2j + infection2_juv, dd_mortality[6])[0];
    double infected2_adult_death = Rcpp::rbinom(1, I2a + infection2_adult, dd_mortality[7])[0];
    double recovery2_juv = std::min(Rcpp::rbinom(1, I2j + infection2_juv - infected2_juvenile_death, recovery[6])[0], I2j + infection2_juv - infected2_juvenile_death);
    double recovery2_adult = std::min(Rcpp::rbinom(1, I2a + infection2_adult - infected2_adult_death, recovery[7])[0], I2a + infection2_adult - infected2_adult_death);
    double recovered_juvenile_death = Rcpp::rbinom(1, Rj + recovery1_juv + recovery2_juv - infection2_juv, dd_mortality[4])[0];
    double recovered_adult_death = Rcpp::rbinom(1, Ra + recovery1_adult + recovery2_adult - infection2_adult, dd_mortality[5])[0];

    // Update state for the next time step
    state(0, t + 1) = isBreedingSeason ? Sj + new_juv - infection1_juv - susceptible_juvenile_death : Sj - infection1_juv - susceptible_juvenile_death;
    state(1, t + 1) = Sa - infection1_adult - susceptible_adult_death;
    state(2, t + 1) = I1j + infection1_juv - recovery1_juv - infected1_juvenile_death;
    state(3, t + 1) = I1a + infection1_adult - recovery1_adult - infected1_adult_death;
    state(4, t + 1) = Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death;
    state(5, t + 1) = Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death;
    state(6, t + 1) = I2j + infection2_juv - recovery2_juv - infected2_juvenile_death;
    state(7, t + 1) = I2a + infection2_adult - recovery2_adult - infected2_adult_death;
  }

  // Return the final state
  return Rcpp::NumericVector(state.col(season_length).begin(), state.col(season_length).end());
}
