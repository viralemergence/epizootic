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
//' @param transmission A vector of length 4 with the transmission rates for each
//' susceptible/recovered stage in the season in question.
//' @param recovery A vector of length 4 with the transmission rates for each
//' susceptible/recovered stage in the season in question.
//' @param fecundity Only necessary when `season = "breeding"` (see below).
//' Default NULL. A single numeric with the daily fecundity of adults.
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
                                   Rcpp::NumericVector recovery, double fecundity,
                             double abundance_threshold, double carrying_capacity,
                             const char * season) {

   int n_stages = 8;
   arma::mat state(n_stages, season_length + 1);
   Rcpp::NumericVector dd_mortality(n_stages);
   const char *str_comp1 = "breeding";

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

     double N = std::min(static_cast<double>(std::accumulate(current_pop_vector.begin(),
                                                             current_pop_vector.end(), 0)),
                                                             carrying_capacity);


     if (N < abundance_threshold) {
       // Fill the rest of the season with 0 and exit
       for (int i = t + 1; i < season_length; i++) {
         for (int j = 0; j < n_stages; j++) {
           state(j, i) = 0.0;
         }
       }
       break;
     }

     dd_mortality = (1.0 + N / carrying_capacity) * mortality;

     for (int i = 0; i < n_stages; i++) {
       dd_mortality[i] = std::min(dd_mortality[i], 1.0);
     }

     double new_juv = 0.0;

     if (std::strcmp(str_comp1, season) == 0) {
       double adults = Sa + I1a + Ra + I2a;
       Rcpp::NumericVector births = Rcpp::rpois(adults, fecundity * (1.0 - N / carrying_capacity));
       new_juv = sum(births);
     }

     Rcpp::NumericVector infection1_juv_raw = Rcpp::rbinom(1, Sj * (I1j + I2j + I1a + I2a), transmission[0]);
     double infection1_juv = std::min(infection1_juv_raw[0], Sj);

     Rcpp::NumericVector infection1_adult_raw = Rcpp::rbinom(1, Sa * (I1j + I2j + I1a + I2a), transmission[1]);
     double infection1_adult = std::min(infection1_adult_raw[0], Sa);

     Rcpp::NumericVector recovery1_juv_raw = Rcpp::rbinom(1, I1j, recovery[0]);
     double recovery1_juv = std::min(recovery1_juv_raw[0], I1j);

     Rcpp::NumericVector recovery1_adult_raw = Rcpp::rbinom(1, I1a, recovery[1]);
     double recovery1_adult = std::min(recovery1_adult_raw[0], I1a);

     Rcpp::NumericVector recovery2_juv_raw = Rcpp::rbinom(1, I2j, recovery[2]);
     double recovery2_juv = std::min(recovery2_juv_raw[0], I2j);

     Rcpp::NumericVector recovery2_adult_raw = Rcpp::rbinom(1, I2a, recovery[3]);
     double recovery2_adult = std::min(recovery2_adult_raw[0], I2a);

     Rcpp::NumericVector infection2_juv_raw = Rcpp::rbinom(1, Rj * (I1j + I2j + I1a + I2a), transmission[2]);
     double infection2_juv = std::min(infection2_juv_raw[0], Rj);

     Rcpp::NumericVector infection2_adult_raw = Rcpp::rbinom(1, Ra * (I1j + I2j + I1a + I2a), transmission[3]);
     double infection2_adult = std::min(infection2_adult_raw[0], Ra);

     Rcpp::NumericVector susceptible_adult_death_raw = Rcpp::rbinom(1, Sa - infection1_adult, dd_mortality[1]);
     double susceptible_adult_death = susceptible_adult_death_raw[0];

     Rcpp::NumericVector recovered_juvenile_death_raw = Rcpp::rbinom(1, Rj + recovery1_juv + recovery2_juv - infection2_juv, dd_mortality[4]);
     double recovered_juvenile_death = recovered_juvenile_death_raw[0];

     Rcpp::NumericVector recovered_adult_death_raw = Rcpp::rbinom(1, Ra + recovery1_adult + recovery2_adult - infection2_adult, dd_mortality[5]);
     double recovered_adult_death = recovered_adult_death_raw[0];

     double susceptible_juvenile_death = 0.0;
     if (std::strcmp(str_comp1, season) == 0) {
       Rcpp::NumericVector susceptible_juvenile_death_raw = Rcpp::rbinom(1, Sj + new_juv - infection1_juv,
                                                     dd_mortality[0]);
       susceptible_juvenile_death = susceptible_juvenile_death_raw[0];
     } else {
       Rcpp::NumericVector susceptible_juvenile_death_raw = Rcpp::rbinom(1, Sj - infection1_juv,
                                                     dd_mortality[0]);
       susceptible_juvenile_death = susceptible_juvenile_death_raw[0];
     }
     Rcpp::NumericVector infected1_juvenile_death_raw = Rcpp::rbinom(1, I1j + infection1_juv - recovery1_juv, dd_mortality[2]);
     double infected1_juvenile_death = infected1_juvenile_death_raw[0];

     Rcpp::NumericVector infected1_adult_death_raw = Rcpp::rbinom(1, I1a + infection1_adult - recovery1_adult, dd_mortality[3]);
     double infected1_adult_death = infected1_adult_death_raw[0];

     Rcpp::NumericVector infected2_juvenile_death_raw = Rcpp::rbinom(1, I2j + infection2_juv - recovery2_juv, dd_mortality[6]);
     double infected2_juvenile_death = infected2_juvenile_death_raw[0];

     Rcpp::NumericVector infected2_adult_death_raw = Rcpp::rbinom(1, I2a + infection2_adult - recovery2_adult, dd_mortality[7]);
     double infected2_adult_death = infected2_adult_death_raw[0];

     // Update state for the next time step
     if (std::strcmp(str_comp1, season) == 0) {
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
