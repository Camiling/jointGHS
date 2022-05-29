#ifndef _jointGHS_ECM_GHS_H
#define _jointGHS_ECM_GHS_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


List ECM_GHS(arma::cube S, arma::cube theta, arma::cube sigma, arma::cube Lambda_sq, arma::vec N, int M, int K, double epsilon, bool verbose, int maxitr, arma::vec tau_sq, bool fix_tau = true)
#endif