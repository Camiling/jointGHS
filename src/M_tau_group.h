#ifndef _fastGHS_M_TAU_GROUP_H
#define _fastGHS_M_TAU_GROUP_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat M_tau_group(int M, mat &theta, mat &Lambda_sq, int exist_group, uvec &group, mat &N_groups, mat E_xi);

#endif