#ifndef _jointGHS_Q_val_H
#define _jointGHS_Q_val_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double Q_val(int N, int M, mat &theta, mat&S, mat &Lambda_sq, mat &E_NuInv,  mat Tau_sq, mat E_XiInv, double tau_sq=0, double E_xiInv=0);

#endif