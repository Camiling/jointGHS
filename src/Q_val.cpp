// [[Rcpp::depends(RcppArmadillo)]]
#include "Q_val.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

double Q_val(int N, int M, mat &theta, mat&S, mat &Lambda_sq, mat &E_NuInv, int exist_group, uvec &group, mat Tau_sq, mat E_XiInv, double tau_sq, double E_xiInv) {
  // If there are groups, Tau_sq and Xi are used. If not, tau_sq and xi are used. If there are no groups, some Tau_sq and Xi must be provided as dummy arguments. 
  double res;
  mat mat_temp;
  mat Tau_G_sq;
  double sum_temp;
  if(exist_group>0){
    Tau_G_sq = Tau_sq.submat(group,group);
    mat_temp = -2*log(Lambda_sq) - pow(theta,2)/Lambda_sq/Tau_G_sq/2 - E_NuInv/sqrt(Lambda_sq) - 2*log(Tau_G_sq) - (1+1/Tau_G_sq)%E_XiInv;
    sum_temp = (sum(sum(mat_temp)) - sum(mat_temp.diag()))/2;
    res = M/2*log(det(theta)) - trace(S*theta)/2 + sum_temp;
  } else{
    mat_temp = -2*log(Lambda_sq) - pow(theta,2)/Lambda_sq/(2*tau_sq) - E_NuInv/sqrt(Lambda_sq);
    sum_temp = (sum(sum(mat_temp)) - sum(mat_temp.diag()))/2;
    res = M/2*log(det(theta)) - trace(S*theta)/2 + sum_temp - (M*(M-1)/2+3)*log(tau_sq)/2 - E_xiInv/tau_sq;
  }
  return res;
}
