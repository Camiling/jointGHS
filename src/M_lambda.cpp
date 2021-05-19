#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>
#include "M_lambda.h"

using namespace std;
using namespace Rcpp;
using namespace arma;


mat M_lambda(int N, int M, mat &theta,mat E_Nu, int exist_group, uvec &group,mat Tau_sq, double tau_sq) {
  // Slighty different computation if there are groups.
  mat Lambda_sq_new;
  if(exist_group>0){
    mat tau_G = Tau_sq.submat(group,group); // MxM mat containing tau of each variable combination
    Lambda_sq_new = (E_Nu + pow(theta,2)/tau_G/2)/2;
  }else{
    Lambda_sq_new = (E_Nu + pow(theta,2)/tau_sq/2)/2;
  }
  return Lambda_sq_new; 
}

