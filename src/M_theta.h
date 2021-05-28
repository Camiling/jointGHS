#ifndef _fastGHS_M_THETA_H
#define _fastGHS_M_THETA_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace arma;
using namespace Rcpp;

cube M_theta(int N, int M, mat theta, mat &S, mat &sigma, mat &Lambda_sq, uvec pseq, int exist_group, uvec &group, mat Tau_sq, double machine_eps, bool stop_underflow, double tau_sq=0);

#endif