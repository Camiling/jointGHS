#ifndef _jointGHS_M_LAMBDA_H
#define _jointGHS_M_LAMBDA_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 
using namespace std;
using namespace Rcpp;
using namespace arma;

mat M_lambda(mat &theta,mat E_Nu, double tau_sq);

#endif