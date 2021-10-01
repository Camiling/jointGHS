#ifndef _jointGHS_E_NU_H
#define _jointGHS_E_NU_H

#include <RcppArmadillo.h>
using namespace arma;
//#include <stdio.h>
//#include <math.h> 
//using namespace std;
//using namespace Rcpp;
//using namespace arma;

mat E_Nu(arma::cube &Lambda_sq, int M, int K);
#endif