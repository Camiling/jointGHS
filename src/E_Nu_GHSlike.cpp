#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

mat E_Nu_GHSlike(double tau_sq, arma::mat theta) {
  // Get the conditional expectation of Nu
  mat ans = (1/log(1+tau_sq/pow(theta, 2))) *tau_sq*tau_sq/((pow(theta,2)+tau_sq) % pow(theta, 2));
  return(ans);
}
