// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// aspatial_siri
Rcpp::NumericVector aspatial_siri(arma::vec initial_pop, int season_length, Rcpp::NumericVector mortality, Rcpp::NumericVector transmission, Rcpp::NumericVector recovery, Rcpp::NumericVector fecundity, double abundance_threshold, double carrying_capacity, const char * season);
RcppExport SEXP _epizootic_aspatial_siri(SEXP initial_popSEXP, SEXP season_lengthSEXP, SEXP mortalitySEXP, SEXP transmissionSEXP, SEXP recoverySEXP, SEXP fecunditySEXP, SEXP abundance_thresholdSEXP, SEXP carrying_capacitySEXP, SEXP seasonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type initial_pop(initial_popSEXP);
    Rcpp::traits::input_parameter< int >::type season_length(season_lengthSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mortality(mortalitySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type transmission(transmissionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type recovery(recoverySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fecundity(fecunditySEXP);
    Rcpp::traits::input_parameter< double >::type abundance_threshold(abundance_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type carrying_capacity(carrying_capacitySEXP);
    Rcpp::traits::input_parameter< const char * >::type season(seasonSEXP);
    rcpp_result_gen = Rcpp::wrap(aspatial_siri(initial_pop, season_length, mortality, transmission, recovery, fecundity, abundance_threshold, carrying_capacity, season));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epizootic_aspatial_siri", (DL_FUNC) &_epizootic_aspatial_siri, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_epizootic(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
