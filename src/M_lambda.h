#ifndef _fastGHS_M_LAMBDA_H
#define _fastGHS_M_LAMBDA_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 
using namespace std;
using namespace Rcpp;
using namespace arma;

mat M_lambda(int N, int M, mat &theta,mat E_Nu, int exist_group, uvec &group,mat Tau_sq, double machine_eps, bool stop_underflow, double tau_sq=0);

#endif