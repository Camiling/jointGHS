#ifndef _fastGHS_E_NU_GHSLIKE_H
#define _fastGHS_E_NU_GHSLIKE_H

#include <RcppArmadillo.h>
using namespace arma;
//#include <stdio.h>
//#include <math.h> 
//using namespace std;
//using namespace Rcpp;
//using namespace arma;

mat E_Nu_GHSlike(double tau_sq, arma::mat theta);
#endif