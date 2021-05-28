#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>
#include "M_lambda.h"

using namespace std;
using namespace Rcpp;
using namespace arma;


mat M_lambda(int N, int M, mat &theta,mat E_Nu, int exist_group, uvec &group,mat Tau_sq, double machine_eps,bool stop_underflow, double tau_sq) {
  // Slighty different computation if there are groups.
  mat Lambda_sq_new;
  if(exist_group>0){
    mat tau_G = Tau_sq.submat(group,group); // MxM mat containing tau of each variable combination
    Lambda_sq_new = (E_Nu + pow(theta,2)/tau_G/2)/2;
  }else{
    Lambda_sq_new = (E_Nu + pow(theta,2)/tau_sq/2)/2;
  }
  int i; 
  int j; 
  
  if(stop_underflow){
    if(min(min(Lambda_sq_new)) < machine_eps){
      for(i=0; i<M; i++){
        for(j=0; j<M; j++){
          if(Lambda_sq_new(i,j) < machine_eps){
            Lambda_sq_new(i,j) = machine_eps;
          }
        }
      }
    } 
  }
  
  return Lambda_sq_new; 
}

