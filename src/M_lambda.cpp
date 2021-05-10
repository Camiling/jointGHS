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
    //mat tau_G = zeros<mat>(M,M); // MxM mat containing tau of each variable combination
    //int i; int j;
    //for(i=0; i < M; i++){
    //  for(j=0; j < M; j++){
    //    tau_G(i,j) = Tau_sq(group(i),group(j));
    //  }
    //}
    mat tau_G = Tau_sq.submat(group,group); // MxM mat containing tau of each variable combination
    Lambda_sq_new = pow((E_Nu + sqrt(pow(E_Nu,2) + 16*pow(theta,2)/tau_G))/8, 2);
  }else{
    Lambda_sq_new = pow((E_Nu + sqrt(pow(E_Nu,2) + 16*pow(theta,2)/tau_sq))/8, 2);
  }
  return Lambda_sq_new; 
}

