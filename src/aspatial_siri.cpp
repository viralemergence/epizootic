#include <RcppArmadillo.h>

//' Helper Function for Seasonal SIRI Simulation
//'
//' Testing out a seasonal version instead of a daily timestep version.
// [[Rcpp::export]]

Rcpp::List aspatial_siri(Rcpp::List pop_list,
                         Rcpp::List mortality_list,
                         Rcpp::List transmission_list,
                         Rcpp::List recovery_list,
                         Rcpp::List fecundity_list,
                         double abundance_threshold,
                         Rcpp::NumericVector carrying_capacity_list,
                         const char * season) {

   int n_stages = 8;
   Rcpp::NumericVector dd_mortality(n_stages);
   Rcpp::NumericVector final_state(n_stages);
   Rcpp::List results(pop_list.size());
   const char *str_comp1 = "breeding";
   Rcpp::NumericVector initial_pop(n_stages);
   double carrying_capacity = 0.0;
   Rcpp::NumericVector mortality(n_stages);
   double fecundity = 0.0;
   Rcpp::NumericVector transmission(n_stages);
   Rcpp::NumericVector recovery(n_stages);
   double susceptible_juvenile_death = 0.0;
   double infection1_juv_raw = 0.0;
   double infection1_adult_raw = 0.0;
   double susceptible_adult_death = 0.0;
   double infected1_juvenile_death = 0.0;
   double infected1_adult_death = 0.0;
   double recovery1_juv_raw = 0.0;
   double recovery1_adult_raw = 0.0;
   double infection2_juv_raw = 0.0;
   double infection2_adult_raw = 0.0;
   double infected2_juvenile_death = 0.0;
   double infected2_adult_death = 0.0;
   double recovery2_juv_raw = 0.0;
   double recovery2_adult_raw = 0.0;
   double recovered_juvenile_death = 0.0;
   double recovered_adult_death = 0.0;

   for (int p = 0; p < pop_list.size(); p++) {
      initial_pop = pop_list[p];
      carrying_capacity = carrying_capacity_list[p];
      mortality = mortality_list[p];
      fecundity = fecundity_list[p];
      transmission = transmission_list[p];
      recovery = recovery_list[p];

      // Unpack initial_pops
      double Sa = initial_pop(1);
      double Sj = initial_pop(0);
      double I1j = initial_pop(2);
      double I1a = initial_pop(3);
      double Rj = initial_pop(4);
      double Ra = initial_pop(5);
      double I2j = initial_pop(6);
      double I2a = initial_pop(7);

      double total_pop = sum(initial_pop);

      double N = std::min(total_pop, carrying_capacity);

      double dd_multiplier = 1.0 + N / carrying_capacity;

      dd_mortality = mortality * std::min(dd_multiplier, 1.0);

      double new_juv = 0.0;

      if (std::strcmp(str_comp1, season) == 0) {
         double adults = Sa + I1a + Ra + I2a;
         Rcpp::NumericVector births = Rcpp::rpois(adults,
                                                  fecundity * (1.0 - N / carrying_capacity));
         new_juv = sum(births);
      }

      infection1_juv_raw = R::rbinom(Sj * (I1j + I2j + I1a + I2a),
                                                            transmission[0]);
      double infection1_juv = std::min(infection1_juv_raw, Sj);

      infection1_adult_raw = R::rbinom(Sa * (I1j + I2j + I1a + I2a),
                                       transmission[1]);
      double infection1_adult = std::min(infection1_adult_raw, Sa);

      susceptible_adult_death = R::rbinom(Sa - infection1_adult, dd_mortality[1]);

      if (std::strcmp(str_comp1, season) == 0) {
         susceptible_juvenile_death = R::rbinom(Sj + new_juv - infection1_juv,
                                                    dd_mortality[0]);
      } else {
         susceptible_juvenile_death = R::rbinom(Sj - infection1_juv,
                                                dd_mortality[0]);
      }

      infected1_juvenile_death = R::rbinom(I1j + infection1_juv,
                                                  dd_mortality[2]);

      infected1_adult_death = R::rbinom(I1a + infection1_adult, dd_mortality[3]);

      recovery1_juv_raw = R::rbinom(I1j + infection1_juv - infected1_juvenile_death,
                                       recovery[0]);
      double recovery1_juv = std::min(recovery1_juv_raw,
                                      I1j + infection1_juv - infected1_juvenile_death);

      recovery1_adult_raw = R::rbinom(I1a + infection1_adult - infected1_adult_death,
                                         recovery[1]);
      double recovery1_adult = std::min(recovery1_adult_raw,
                                        I1a + infection1_adult - infected1_adult_death);

      infection2_juv_raw = R::rbinom((Rj + recovery1_juv) * (I1j + I2j + I1a + I2a),
                                                            transmission[2]);
      double infection2_juv = std::min(infection2_juv_raw, Rj + infection1_juv);

      infection2_adult_raw = R::rbinom((Ra + recovery1_adult) * (I1j + I2j + I1a + I2a),
                                                              transmission[3]);
      double infection2_adult = std::min(infection2_adult_raw,
                                         Ra + recovery1_adult);

      infected2_juvenile_death = R::rbinom(I2j + infection2_juv,
                                               dd_mortality[6]);

      infected2_adult_death = R::rbinom(I2a + infection2_adult,
                                            dd_mortality[7]);

      recovery2_juv_raw = R::rbinom(I2j + infection2_juv - infected2_juvenile_death,
                                                           recovery[2]);
      double recovery2_juv = std::min(recovery2_juv_raw, I2j + infection2_juv - infected2_juvenile_death);

      recovery2_adult_raw = R::rbinom(I2a + infection2_adult - infected2_adult_death,
                                                             recovery[3]);
      double recovery2_adult = std::min(recovery2_adult_raw,
                                        I2a + infection2_adult - infected2_juvenile_death);

      recovered_juvenile_death = R::rbinom(Rj + recovery1_juv + recovery2_juv - infection2_juv,
                                           dd_mortality[4]);

      recovered_adult_death = R::rbinom(Ra + recovery1_adult + recovery2_adult - infection2_adult,
                                            dd_mortality[5]);

      // Update state for the next time step
      if (std::strcmp(str_comp1, season) == 0) {
         final_state(0) = Sj + new_juv - infection1_juv - susceptible_juvenile_death;
      } else {
         final_state(0) = Sj - infection1_juv - susceptible_juvenile_death;
      }
      final_state(1) = Sa - infection1_adult - susceptible_adult_death;
      final_state(2) = I1j + infection1_juv - recovery1_juv -
         infected1_juvenile_death;
      final_state(3) = I1a + infection1_adult  - recovery1_adult -
         infected1_adult_death;
      final_state(4) = Rj + recovery1_juv + recovery2_juv - infection2_juv -
         recovered_juvenile_death;
      final_state(5) = Ra + recovery1_adult + recovery2_adult -
         infection2_adult - recovered_adult_death;
      final_state(6) = I2j + infection2_juv - recovery2_juv -
         infected2_juvenile_death;
      final_state(7) = I2a + infection2_adult - recovery2_adult -
         infected2_adult_death;

      if (N < abundance_threshold) {
         // Fill everything with 0 and exit
         for (int j = 0; j < n_stages; j++) {
            final_state(j) = 0.0;
         }
      }

      results[p] = final_state;
   }

   return results;
}
