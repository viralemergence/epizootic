#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Helper Function for Seasonal SIR Simulation
//'
//' This function is an internal one that does the aspatial simulations within
//' one population for one timestep, for any given season. It is the two-stage
//' (juvenile, adult) SIR analogue of [aspatial_siri()], in which recovery
//' confers permanent immunity.
//'
//' @name aspatial_sir
//'
//' @param initial_pop A vector of length 6 showing the initial abundance for
//' each combination of stage and compartment, in the order susceptible
//' juvenile, susceptible adult, infected juvenile, infected adult, recovered
//' juvenile, recovered adult.
//' @param season_length The length of the season in days.
//' @param mortality A vector of length 6 with the mortality rates for each
//' stage and compartment in the season in question.
//' @param transmission A vector of length 6 with the transmission rates for
//' each stage and compartment in the season in question. Only the susceptible
//' entries (elements 1 and 2) are used.
//' @param recovery A vector of length 6 with the recovery rates for each
//' stage and compartment in the season in question. Only the infected entries
//' (elements 3 and 4) are used.
//' @param fecundity A vector of length 6 with the fecundity for each
//' reproductive segment. Only adult fecundity (element 2) is used.
//' @param abundance_threshold A quasi-extinction threshold below which a
//' population becomes extinct.
//' @param carrying_capacity A single numeric that indicates the carrying
//' capacity of the population in this season.
//' @param season Either "breeding" or "non-breeding."
//' @return A vector of length 6 showing the abundance for each combination of
//' stage and compartment at the end of the season.
//'
//' @examples
//' aspatial_sir(
//'  initial_pop = c(50000, 50000, 1, 0, 0, 0),
//'  season_length = 100,
//'  mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0),
//'  fecundity = c(0, 15/182, 0, 15/182, 0, 15/182),
//'  transmission = c(0.00002, 0.00001, 0, 0, 0, 0),
//'  recovery = c(0, 0, 0.05714286, 0.05714286, 0, 0),
//'  carrying_capacity = 150000,
//'  abundance_threshold = 10,
//'  season = "breeding"
//' )
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector aspatial_sir(Rcpp::NumericVector initial_pop,
                                 int season_length,
                                 Rcpp::NumericVector mortality,
                                 Rcpp::NumericVector transmission,
                                 Rcpp::NumericVector recovery,
                                 Rcpp::NumericVector fecundity,
                                 double abundance_threshold,
                                 double carrying_capacity,
                                 const std::string& season) {

  int n_stages = 6;
  arma::mat state(n_stages, season_length + 1);
  Rcpp::NumericVector dd_mortality(n_stages);
  bool isBreedingSeason = (season == "breeding");
  double birth_rate = fecundity[1];

  state.col(0) = Rcpp::as<arma::vec>(initial_pop);
  for (int t = 0; t < season_length; t++) {
    // Unpack states
    double Sj = state(0, t);
    double Sa = state(1, t);
    double Ij = state(2, t);
    double Ia = state(3, t);
    double Rj = state(4, t);
    double Ra = state(5, t);

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
      double adults = Sa + Ia + Ra;
      double birth_rate_adj = birth_rate * (1.0 - N / carrying_capacity);
      new_juv = Rcpp::rpois(1, adults * birth_rate_adj)[0];
    }

    bool hasInfections = (Ij + Ia) > 0;

    if (!hasInfections) {
      // Only update susceptibles and recovereds, skipping all infection and
      // recovery calculations
      state(0, t + 1) = isBreedingSeason ? Sj + new_juv - Rcpp::rbinom(1, Sj + new_juv, dd_mortality[0])[0] : Sj - Rcpp::rbinom(1, Sj, dd_mortality[0])[0];
      state(1, t + 1) = Sa - Rcpp::rbinom(1, Sa, dd_mortality[1])[0];
      state(2, t + 1) = Ij;
      state(3, t + 1) = Ia;
      state(4, t + 1) = Rj - Rcpp::rbinom(1, Rj, dd_mortality[4])[0];
      state(5, t + 1) = Ra - Rcpp::rbinom(1, Ra, dd_mortality[5])[0];
      continue;
    }

    double infection_juv = std::min(Rcpp::rbinom(1, Sj * (Ij + Ia), transmission[0])[0], Sj);
    double infection_adult = std::min(Rcpp::rbinom(1, Sa * (Ij + Ia), transmission[1])[0], Sa);
    double susceptible_adult_death = Rcpp::rbinom(1, Sa - infection_adult, dd_mortality[1])[0];
    double susceptible_juvenile_death = Rcpp::rbinom(1, isBreedingSeason ? Sj + new_juv - infection_juv : Sj - infection_juv, dd_mortality[0])[0];
    double infected_juvenile_death = Rcpp::rbinom(1, Ij + infection_juv, dd_mortality[2])[0];
    double infected_adult_death = Rcpp::rbinom(1, Ia + infection_adult, dd_mortality[3])[0];
    double recovery_juv = std::min(Rcpp::rbinom(1, Ij + infection_juv - infected_juvenile_death, recovery[2])[0], Ij + infection_juv - infected_juvenile_death);
    double recovery_adult = std::min(Rcpp::rbinom(1, Ia + infection_adult - infected_adult_death, recovery[3])[0], Ia + infection_adult - infected_adult_death);
    double recovered_juvenile_death = Rcpp::rbinom(1, Rj + recovery_juv, dd_mortality[4])[0];
    double recovered_adult_death = Rcpp::rbinom(1, Ra + recovery_adult, dd_mortality[5])[0];

    // Update state for the next time step
    state(0, t + 1) = isBreedingSeason ? Sj + new_juv - infection_juv - susceptible_juvenile_death : Sj - infection_juv - susceptible_juvenile_death;
    state(1, t + 1) = Sa - infection_adult - susceptible_adult_death;
    state(2, t + 1) = Ij + infection_juv - recovery_juv - infected_juvenile_death;
    state(3, t + 1) = Ia + infection_adult - recovery_adult - infected_adult_death;
    state(4, t + 1) = Rj + recovery_juv - recovered_juvenile_death;
    state(5, t + 1) = Ra + recovery_adult - recovered_adult_death;
  }

  // Return the final state
  return Rcpp::NumericVector(state.col(season_length).begin(), state.col(season_length).end());
}
