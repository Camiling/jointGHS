#include "E_xi.h"

double E_xi(double &tau_sq) {
  // Get the conditional expectation of 1/xi. 
  // No grouping
  double ans = tau_sq/(tau_sq+1);
  return ans;
}