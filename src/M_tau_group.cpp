// [[Rcpp::depends(RcppArmadillo)]]
#include "M_tau_group.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

mat M_tau_group(int M, mat &theta, mat &Lambda_sq, int exist_group, uvec &group, mat &N_groups, mat E_xi) {
  // maximize Tau_sq
  // With grouping, so E_xi is now exist_group x exist_group
  mat mat_temp = pow(theta,2)/Lambda_sq;
  mat sum_temp = zeros<mat>(exist_group,exist_group); // ngroups x ngroups 
  int i;
  int j;
  for(i=0; i < M; i++){
    for(j=0; j < M; j++){
      if(i < j){ 
        sum_temp(group(i), group(j)) += mat_temp(i,j);
        if(group(i)!=group(j)){ // If not same group, add to other side of the diagonal too for symmetry. 
          sum_temp(group(j), group(i)) += mat_temp(i,j);
        }
      }
    }
  }
  // N_groups is an ngroups x ngroups mat with the number of variables in each combination group
  mat tau_sq = (sum_temp + 2*E_xi%N_groups)/(4*N_groups);
  return tau_sq; 
}