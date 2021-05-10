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
  if(exist_group>0){
    
  }else{
    
  }
  return 1;
}
