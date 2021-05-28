// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ECM_GHS
List ECM_GHS(arma::mat X, arma::mat S, arma::mat theta, arma::mat sigma, arma::mat Lambda_sq, double epsilon, bool verbose, int maxitr, bool savepath, int exist_group, arma::uvec group, arma::mat N_groups, bool save_Q, double tau_sq, arma::mat Tau_sq, double machine_eps, bool use_ICM, bool stop_underflow);
RcppExport SEXP _fastGHS_ECM_GHS(SEXP XSEXP, SEXP SSEXP, SEXP thetaSEXP, SEXP sigmaSEXP, SEXP Lambda_sqSEXP, SEXP epsilonSEXP, SEXP verboseSEXP, SEXP maxitrSEXP, SEXP savepathSEXP, SEXP exist_groupSEXP, SEXP groupSEXP, SEXP N_groupsSEXP, SEXP save_QSEXP, SEXP tau_sqSEXP, SEXP Tau_sqSEXP, SEXP machine_epsSEXP, SEXP use_ICMSEXP, SEXP stop_underflowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda_sq(Lambda_sqSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    Rcpp::traits::input_parameter< bool >::type savepath(savepathSEXP);
    Rcpp::traits::input_parameter< int >::type exist_group(exist_groupSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type N_groups(N_groupsSEXP);
    Rcpp::traits::input_parameter< bool >::type save_Q(save_QSEXP);
    Rcpp::traits::input_parameter< double >::type tau_sq(tau_sqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Tau_sq(Tau_sqSEXP);
    Rcpp::traits::input_parameter< double >::type machine_eps(machine_epsSEXP);
    Rcpp::traits::input_parameter< bool >::type use_ICM(use_ICMSEXP);
    Rcpp::traits::input_parameter< bool >::type stop_underflow(stop_underflowSEXP);
    rcpp_result_gen = Rcpp::wrap(ECM_GHS(X, S, theta, sigma, Lambda_sq, epsilon, verbose, maxitr, savepath, exist_group, group, N_groups, save_Q, tau_sq, Tau_sq, machine_eps, use_ICM, stop_underflow));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastGHS_ECM_GHS", (DL_FUNC) &_fastGHS_ECM_GHS, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastGHS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
