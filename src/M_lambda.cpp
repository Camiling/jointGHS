// [[Rcpp::depends(RcppArmadillo)]]
#include "M_lambda.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;


mat M_lambda(mat &theta, mat E_Nu, double tau_sq) {
  mat Lambda_sq_new;
  Lambda_sq_new = (E_Nu + pow(theta,2)/tau_sq/2)/2;
  return Lambda_sq_new; 
}

