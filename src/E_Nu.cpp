#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>
//#include "E_Nu.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

mat E_Nu(List &Lambda_sq) {
  // Get the conditional expectation of 1/Nu
  mat ans = Lambda_sq/(Lambda_sq+1);
  return(ans);
}
