// [[Rcpp::depends(RcppArmadillo)]]
#include "M_tau.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

double M_tau(int M, mat &theta, mat &Lambda_sq, double E_xi) {
  // maximize Tau_sq
  // No grouping
  // Maybe Rewrite: sum theta_ij/lambda_ij to sum exp (log theta_ij^2 - log lambda_ij^2)
  // double sum_temp = (exp(log(pow(theta,2))-log(Lambda_sq)) - );
  mat mat_temp = pow(theta,2)/Lambda_sq;
  double sum_temp = (sum(sum(mat_temp)) - sum(mat_temp.diag()))/2;
  double tau_sq = (2*sum_temp + 4*E_xi)/(M*(M-1)+6);
  return tau_sq;
}