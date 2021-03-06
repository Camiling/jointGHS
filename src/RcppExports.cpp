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
List ECM_GHS(arma::cube S, arma::cube theta, arma::cube sigma, arma::cube Lambda_sq, arma::vec N, int M, int K, double epsilon, bool verbose, int maxitr, arma::vec tau_sq, bool fix_tau);
RcppExport SEXP _jointGHS_ECM_GHS(SEXP SSEXP, SEXP thetaSEXP, SEXP sigmaSEXP, SEXP Lambda_sqSEXP, SEXP NSEXP, SEXP MSEXP, SEXP KSEXP, SEXP epsilonSEXP, SEXP verboseSEXP, SEXP maxitrSEXP, SEXP tau_sqSEXP, SEXP fix_tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Lambda_sq(Lambda_sqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau_sq(tau_sqSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_tau(fix_tauSEXP);
    rcpp_result_gen = Rcpp::wrap(ECM_GHS(S, theta, sigma, Lambda_sq, N, M, K, epsilon, verbose, maxitr, tau_sq, fix_tau));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_jointGHS_ECM_GHS", (DL_FUNC) &_jointGHS_ECM_GHS, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_jointGHS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
