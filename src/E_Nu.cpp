#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>
//#include "E_Nu.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

mat E_Nu(arma::cube &Lambda_sq, int M, int K) {
  // Get the conditional expectation of 1/Nu
  cube lam_prod = 1/Lambda_sq;
  mat lam_sum = zeros<mat>(M,M);
  int k;
  for (k=0; k < K; k++){
    lam_sum = lam_sum + lam_prod.slice(k);
  }
  mat ans = K/2/(1+lam_sum);
  return(ans);
}
