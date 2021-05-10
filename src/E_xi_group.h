#ifndef _fastGHS_E_XI_GROUP_H
#define _fastGHS_E_XI_GROUP_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace arma;
using namespace Rcpp;

mat E_xi_group(mat &tau_sq);

#endif