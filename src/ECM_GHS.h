#ifndef _jointGHS_ECM_GHS_H
#define _jointGHS_ECM_GHS_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


List ECM_GHS(List X, List S, List theta, List sigma, List Lambda_sq, arma::vec N, int M, int K, double epsilon, bool verbose, int maxitr, bool savepath, bool save_Q, arma::vec tau_sq, bool use_ICM = false, bool fix_tau = false)
#endif