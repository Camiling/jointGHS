#include "E_xi_group.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

mat E_xi_group(mat &tau_sq) {
  // Get the conditional expectations of 1/xi_{gi, gj}
  // Grouped variables 
  // tau_sq is ngroup x ngroup 
  mat ans = tau_sq/(tau_sq+1);
  return(ans);
}