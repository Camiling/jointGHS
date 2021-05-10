#ifndef _fastGHS_Q_val_H
#define _fastGHS_Q_val_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double Q_val(int N, int M, mat &theta, mat&S, mat &Lambda_sq, mat &E_NuInv, int exist_group, uvec &group, mat Tau_sq, mat E_XiInv, double tau_sq=0, double E_xiInv=0);

#endif