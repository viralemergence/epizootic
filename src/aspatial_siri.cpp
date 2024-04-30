#include <RcppArmadillo.h>
#include <random>
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
//' @export
// [[Rcpp::export]]

 Rcpp::NumericVector aspatial_siri(arma::vec initial_pop,
                                   int season_length,
                                   Rcpp::NumericVector mortality,
                                   Rcpp::NumericVector transmission,
                                   Rcpp::NumericVector recovery,
                                   Rcpp::NumericVector fecundity,
                                   double abundance_threshold,
                                   double carrying_capacity,
                                   const char * season) {
   std::random_device rd;
   std::mt19937 gen(rd());
   int n_stages = 8;
   arma::mat state(n_stages, season_length + 1);
   Rcpp::NumericVector dd_mortality(n_stages);
   const char *str_comp1 = "breeding";
   bool isBreedingSeason = std::strcmp(str_comp1, season) == 0;
   double birth_rate = fecundity[1];

   state.col(0) = initial_pop;
   for (int t = 0; t < season_length; t++) {
     // Unpack states
     double Sa = state(1, t);
     double Sj = state(0, t);
     double I1j = state(2, t);
     double I1a = state(3, t);
     double Rj = state(4, t);
     double Ra = state(5, t);
     double I2j = state(6, t);
     double I2a = state(7, t);

     arma::subview_col<double> current_pop = state.col(t);
     std::vector<double> current_pop_vector(current_pop.begin(), current_pop.end());

     double N = std::min(arma::accu(state.col(t)), carrying_capacity);

     if (N < abundance_threshold) {
      state.cols(t + 1, season_length).zeros();
      break;
     }

     dd_mortality = (1.0 + N / carrying_capacity) * mortality;

     for (int i = 0; i < n_stages; i++) {
       dd_mortality[i] = std::min(dd_mortality[i], 1.0);
     }

     int new_juv = 0;

     if (isBreedingSeason) {
      double adults = Sa + I1a + Ra + I2a;
      std::poisson_distribution<int> poisson_dist(birth_rate * (1.0 - N / carrying_capacity));

      int total_births = 0;
      for (int i = 0; i < static_cast<int>(adults); ++i) {
          total_births += poisson_dist(gen);
      }

      new_juv = total_births;
    }

    int trials_infection1_juv = static_cast<int>(Sj * (I1j + I2j + I1a + I2a));
    std::binomial_distribution<int> dist_infection1_juv(trials_infection1_juv, transmission[0]);
    double infection1_juv = std::min(static_cast<double>(dist_infection1_juv(gen)), Sj);

    int trials_infection1_adult = static_cast<int>(Sa * (I1j + I2j + I1a + I2a));
    std::binomial_distribution<int> dist_infection1_adult(trials_infection1_adult, transmission[1]);
    double infection1_adult = std::min(static_cast<double>(dist_infection1_adult(gen)), Sa);

    int trials_susceptible_adult_death = static_cast<int>(Sa - infection1_adult);
    std::binomial_distribution<int> dist_susceptible_adult_death(trials_susceptible_adult_death, dd_mortality[1]);
    double susceptible_adult_death = dist_susceptible_adult_death(gen);

    int trials_susceptible_juvenile_death = isBreedingSeason ? static_cast<int>(Sj + new_juv - infection1_juv) : static_cast<int>(Sj - infection1_juv);
    std::binomial_distribution<int> dist_susceptible_juvenile_death(trials_susceptible_juvenile_death, dd_mortality[0]);
    double susceptible_juvenile_death = dist_susceptible_juvenile_death(gen);

    // Infected juvenile deaths
    int trials_infected1_juvenile = static_cast<int>(I1j + infection1_juv);
    std::binomial_distribution<int> dist_infected1_juvenile(trials_infected1_juvenile, dd_mortality[2]);
    double infected1_juvenile_death = dist_infected1_juvenile(gen);

    // Infected adult deaths
    int trials_infected1_adult = static_cast<int>(I1a + infection1_adult);
    std::binomial_distribution<int> dist_infected1_adult(trials_infected1_adult, dd_mortality[3]);
    double infected1_adult_death = dist_infected1_adult(gen);

    // Recovery of juvenile after first infection
    int trials_recovery1_juv = static_cast<int>(I1j + infection1_juv - infected1_juvenile_death);
    std::binomial_distribution<int> dist_recovery1_juv(trials_recovery1_juv, recovery[2]);
    double recovery1_juv = std::min(static_cast<double>(dist_recovery1_juv(gen)),
                                    I1j + infection1_juv - infected1_juvenile_death);

    // Recovery of adult after first infection
    int trials_recovery1_adult = static_cast<int>(I1a + infection1_adult - infected1_adult_death);
    std::binomial_distribution<int> dist_recovery1_adult(trials_recovery1_adult, recovery[3]);
    double recovery1_adult = std::min(static_cast<double>(dist_recovery1_adult(gen)),
                                      I1a + infection1_adult - infected1_adult_death);

    // Second infection juvenile
    int trials_infection2_juv = static_cast<int>((Rj + recovery1_juv) * (I1j + I2j + I1a + I2a));
    std::binomial_distribution<int> dist_infection2_juv(trials_infection2_juv, transmission[4]);
    double infection2_juv = std::min(static_cast<double>(dist_infection2_juv(gen)),
                                     Rj + recovery1_juv);

    // Second infection adult
    int trials_infection2_adult = static_cast<int>((Ra + recovery1_adult) * (I1j + I2j + I1a + I2a));
    std::binomial_distribution<int> dist_infection2_adult(trials_infection2_adult, transmission[5]);
    double infection2_adult = std::min(static_cast<double>(dist_infection2_adult(gen)),
                                       Ra + recovery1_adult);

    // Second infected juvenile deaths
    int trials_infected2_juvenile = static_cast<int>(I2j + infection2_juv);
    std::binomial_distribution<int> dist_infected2_juvenile(trials_infected2_juvenile, dd_mortality[6]);
    double infected2_juvenile_death = dist_infected2_juvenile(gen);

    // Second infected adult deaths
    int trials_infected2_adult = static_cast<int>(I2a + infection2_adult);
    std::binomial_distribution<int> dist_infected2_adult(trials_infected2_adult, dd_mortality[7]);
    double infected2_adult_death = dist_infected2_adult(gen);

    // Second recovery juvenile
    int trials_recovery2_juv = static_cast<int>(I2j + infection2_juv - infected2_juvenile_death);
    std::binomial_distribution<int> dist_recovery2_juv(trials_recovery2_juv, recovery[6]);
    double recovery2_juv = std::min(static_cast<double>(dist_recovery2_juv(gen)), I2j + infection2_juv - infected2_juvenile_death);

    // Second recovery adult
    int trials_recovery2_adult = static_cast<int>(I2a + infection2_adult - infected2_adult_death);
    std::binomial_distribution<int> dist_recovery2_adult(trials_recovery2_adult, recovery[7]);
    double recovery2_adult = std::min(static_cast<double>(dist_recovery2_adult(gen)), I2a + infection2_adult - infected2_adult_death);

    // Recovered juvenile deaths
    int trials_recovered_juvenile_death = static_cast<int>(Rj + recovery1_juv + recovery2_juv - infection2_juv);
    std::binomial_distribution<int> dist_recovered_juvenile_death(trials_recovered_juvenile_death, dd_mortality[4]);
    double recovered_juvenile_death = dist_recovered_juvenile_death(gen);

    // Recovered adult deaths
    int trials_recovered_adult_death = static_cast<int>(Ra + recovery1_adult + recovery2_adult - infection2_adult);
    std::binomial_distribution<int> dist_recovered_adult_death(trials_recovered_adult_death, dd_mortality[5]);
    double recovered_adult_death = dist_recovered_adult_death(gen);

     // Update state for the next time step
     if (isBreedingSeason) {
       state(0, t + 1) = Sj + new_juv - infection1_juv - susceptible_juvenile_death;
     } else {
       state(0, t + 1) = Sj - infection1_juv - susceptible_juvenile_death;
     }
     state(1, t + 1) = Sa - infection1_adult - susceptible_adult_death;
     state(2, t + 1) = I1j + infection1_juv - recovery1_juv -
       infected1_juvenile_death;
     state(3, t + 1) = I1a + infection1_adult  - recovery1_adult -
         infected1_adult_death;
     state(4, t + 1) = Rj + recovery1_juv + recovery2_juv - infection2_juv -
           recovered_juvenile_death;
     state(5, t + 1) = Ra + recovery1_adult + recovery2_adult -
             infection2_adult - recovered_adult_death;
     state(6, t + 1) = I2j + infection2_juv - recovery2_juv -
       infected2_juvenile_death;
     state(7, t + 1) = I2a + infection2_adult - recovery2_adult -
       infected2_adult_death;
   }

   // Return the final state
   Rcpp::NumericVector final_state(n_stages);
   for (int i = 0; i < n_stages; i++) {
     final_state[i] = state(i, season_length);
   }

   return final_state;
 }
