#ifndef _fastGHS_M_TAU_H
#define _fastGHS_M_TAU_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double M_tau(int M, mat &theta, mat &Lambda_sq, double E_xi, double machine_eps, bool stop_underflow) ;

#endif