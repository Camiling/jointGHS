#ifndef _jointGHS_M_THETA_H
#define _jointGHS_M_THETA_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace arma;
using namespace Rcpp;

cube M_theta(int N, int M, mat theta, mat S, mat sigma, mat Lambda_sq, uvec pseq, double tau_sq);

#endif