#ifndef _jointGHS_ECM_GHS_H
#define _jointGHS_ECM_GHS_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


List ECM_GHS(arma:: S, arma:: theta, arma:: sigma, arma:: Lambda_sq, arma::vec N, int M, int K, double epsilon, bool verbose, int maxitr, bool savepath, bool save_Q, arma::vec tau_sq, bool use_ICM = false, bool fix_tau = true)
#endif