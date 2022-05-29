//#include "ECM_GHS.h"
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

#include "E_xi.h"
#include "M_theta.h"
#include "M_tau.h"
#include "E_Nu.h"
#include "M_lambda.h"

using namespace std;
using namespace Rcpp;
using namespace arma;



//' Perform jointGHS with ECM
//' 
//' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the joint graphical horseshoe
//' 
//' @param S cube of the \eqn{p} by \eqn{p} scatter matrices of the \eqn{K} data sets
//' @param theta cube of the \eqn{K} initial \eqn{p} by \eqn{p} precision matrices
//' @param sigma cube of the \eqn{K} initial \eqn{p} by \eqn{p} covariance matrices
//' @param Lambda_sq cube of the \eqn{K} initial \eqn{p} by \eqn{p} matrices of squared local shrinkage parameters
//' @param N vector of length \eqn{K} of the number of obervations in each data set. 
//' @param M the number of variables/nodes. Integer
//' @param K the number of networks. Integer
//' @param epsilon tolerance for the convergence assessment
//' @param verbose logical indicator of printing information at each iteration
//' @param maxitr maximum number of iterations
//' @param tau_sq vector of initial values of the \eqn{K} squared global shrinkage parameters
//' @param fix_tau logical. Should tau be fixed?
//' 
//' @return A List with resulting ECM estimates, and saved path and objective function convergence information if requested
//' 
//' @export
//'
// [[Rcpp::export]]
List ECM_GHS(arma::cube S, arma::cube theta, arma::cube sigma, arma::cube Lambda_sq, arma::vec N, int M, int K, double epsilon, bool verbose, int maxitr, arma::vec tau_sq, bool fix_tau = true) {


  // initialize intermediate values
  int niter,i, count, k;
  arma::vec eps(K);
  arma::uvec pseq(M);
  for(i = 0; i < M; i++){
    pseq(i) = i;
  }
  
  // Initialise updates
 arma::cube theta_update = theta; // Save previous estimate to assess convergence
 arma::mat E_NuInv(M,M);
 double E_xiInv; // Reused for all graphs
 arma::cube theta_sigma_update(M,M,2);
 eps = eps + epsilon + 1; // Make the initial value sufficiently large
 niter = 1;
 count = 0;
 List list;
 
  while(max(eps)>epsilon & count < maxitr){
      
    // E-step
    E_NuInv = E_Nu(Lambda_sq, M, K); // Found using all K graphs
    
    // For each network separately
    for(k=0; k < K; k++){
      // E-step
      E_xiInv = E_xi(tau_sq(k)); // Only used below, so can use the same variable for all K graphs
      
      // CM-step
      // For standard GHS, we update Lambda_sq, tau and theta in the CM-step
      Lambda_sq.slice(k) = M_lambda(theta.slice(k), E_NuInv, tau_sq(k));
      if (fix_tau == false){
        tau_sq(k) = M_tau(M, theta.slice(k), Lambda_sq.slice(k), E_xiInv);
      }
      theta_sigma_update = M_theta(N(k), M, theta.slice(k), S.slice(k), sigma.slice(k), Lambda_sq.slice(k), pseq, tau_sq(k)); // Pass S as dummy argument
      theta_update.slice(k) = theta_sigma_update.slice(0);
      sigma.slice(k) = theta_sigma_update.slice(1);
      
      eps(k) = max(max(abs(theta_update.slice(k) - theta.slice(k))));
      
      theta.slice(k) = theta_update.slice(k);
    }
    
    count++;
    if(verbose){
      Rcout << "Itr = " << count << " Max diff = " << max(eps) << endl;
    }else{
      Rcout << ".";
    }
  }
  theta = theta_update; // update whole list
  Rcout << " done" << endl;

  // Save results
  list["S"] = S;
  list["theta"] = theta;  
  list["sigma"] = sigma;
  list["Lambda_sq"] = Lambda_sq;
  list["E_NuInv"] = E_NuInv;
  list["tau_sq"] = tau_sq;
  
  return list;
}