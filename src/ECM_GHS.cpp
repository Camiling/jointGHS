//#include "ECM_GHS.h"
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

#include "E_xi.h"
#include "E_xi_group.h"
#include "M_theta.h"
#include "M_tau.h"
#include "M_tau_group.h"
#include "E_Nu.h"
#include "M_lambda.h"

#include "Q_val.h"

using namespace std;
using namespace Rcpp;
using namespace arma;




//' Perform GHS with ECM
//' 
//' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the graphical horseshoe
//' 
//' @param X \eqn{n} by \eqn{p} matrix of data
//' @param S \eqn{p} by \eqn{p} scatter matrix of the data
//' @param theta initial value of the \eqn{p} by \eqn{p} precision matrix
//' @param sigma initial value of the \eqn{p} by \eqn{p} covariance matrix
//' @param Lambda_sq initial value of matrix of squared local shrinkage parameters
//' @param epsilon tolerance for the convergence assessment
//' @param verbose logical indicator of printing information at each iteration
//' @param maxitr maximum number of iterations
//' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm
//' @param exist_group logical. Are the variables grouped?
//' @param group grouping information.
//' @param N_groups If exist_group==T, the number of groups
//' @param save_Q should the value of the objective function at each step be saved?
//' @param tau_sq initial value of squared global shrinkage parameter. If exist_group==T, a dummy value should be provided
//' @param Tau_sq if exist_group==T, an \eqn{ngroup} by \eqn{ngroup} matrix of initial values of the squared global shrinkage parameters within and between groups. If exist_group==F, a dummy value should be provided
//' 
//' @return A List with resulting ECM estimates, and saved path and objective function convergence information if requested
//' 
//' @export
// [[Rcpp::export]]
List ECM_GHS(arma::mat X, arma::mat S, arma::mat theta, arma::mat sigma, arma::mat Lambda_sq, double epsilon, bool verbose, int maxitr, bool savepath, int exist_group, arma::uvec group, arma::mat N_groups, bool save_Q, double tau_sq, arma::mat Tau_sq) {
  // organize input
  //mat X = as<mat>(X_r);
  //mat S = as<mat>(S_r);
  //mat theta = as<mat>(theta_r);
  //mat sigma= as<mat>(sigma_r);
  //mat Lambda_sq = as<mat>(Lambda_sq_r);
  //int exist_group = as<int>(exist_group_r);
  //mat Tau_sq = as<mat>(tau_sq_r); // only used if variables are grouped
  //double tau_sq = as<double>(tau_sq_r); // Ungrouped variables
  
  //double epsilon = as<double>(epsilon_r);
  //int maxitr = as<int>(maxitr_r);
  //bool verbose = as<bool>(verbose_r);
  //bool savepath = as<bool>(savepath_r);
  //bool save_Q = as<bool>(save_Q_r);
  //mat N_groups = as<mat>(N_groups_r);
  //uvec group = as<uvec>(group_r);;
  
  
  // get dimensions
  const int M = X.n_cols;
  const int N = X.n_rows;
  
  // For saving variables
  arma::cube theta_path(M, M, maxitr);
  arma::uvec Q_vals(maxitr);
  
  // initialize intermediate values
  int niter,i, count;
  double eps;
  arma::uvec pseq(M);
  for(i = 0; i < M; i++){
    pseq(i) = i;
  }
  
  // Initialise updates
 arma::mat theta_update = theta; // Save previous estimate to assess convergence
 arma::mat E_NuInv(M,M);
 arma::mat E_XiInv(M,M); // Used if variables are grouped
  double E_xiInv; // Used if variables are not grouped
  arma::cube theta_sigma_update(M,M,2);
  eps = epsilon + 1;
  niter = 1;
  count = 0;
  
  if(savepath){
    theta_path.slice(0) = theta;
  }
    
    
  while(eps>epsilon & count < maxitr){
      
    // E-step
    E_NuInv = E_Nu(Lambda_sq);
    if (exist_group>0){
      E_XiInv = E_xi_group(Tau_sq);
    }
    else {
      E_xiInv = E_xi(tau_sq);
    }
      
    // M-step
    if (exist_group>0){
      Lambda_sq = M_lambda(N, M, theta, E_NuInv, exist_group, group, Tau_sq);
      Tau_sq = M_tau_group(M, theta, Lambda_sq, exist_group, group, N_groups, E_XiInv);
      theta_sigma_update = M_theta(N, M, theta, S, sigma, Lambda_sq, pseq, exist_group, group, Tau_sq);
    }
    else{
      Lambda_sq = M_lambda(N, M, theta, E_NuInv, exist_group, group, S, tau_sq);
      tau_sq = M_tau(M, theta, Lambda_sq, E_xiInv);
      theta_sigma_update = M_theta(N, M, theta, S, sigma, Lambda_sq, pseq, exist_group, group, S, tau_sq); // Pass S as dummy argument
    }
    theta_update = theta_sigma_update.slice(0);
    sigma = theta_sigma_update.slice(1);
    
    if(savepath){
      theta_path.slice(count+1) = theta_update;
    }
    if(save_Q){
      if(exist_group>0){
        Q_vals(count) = Q_val(N, M, theta, S, Lambda_sq, E_NuInv, exist_group, group, Tau_sq, E_XiInv); 
      }else {
        Q_vals(count) = Q_val(N, M, theta, S, Lambda_sq, E_NuInv, exist_group, group, S, S, tau_sq, E_xiInv);  // Pass S as dummy argument
      }
    }
    
    eps = max(max(abs(theta_update - theta)));
    theta = theta_update;
    count++;
    if(verbose){
      Rcout << "Itr = " << count << " Max diff = " << eps << endl;
    }else{
      Rcout << ".";
    }
  }
  theta = theta_update;
  Rcout << " done" << endl;
  
  // Save results
  List list;
  list["S"] = S;
  list["theta"] = theta;  
  list["sigma"] = sigma;
  list["Lambda_sq"] = Lambda_sq;
  if(exist_group>0){
    list["tau_sq"] = Tau_sq;
  }
  else{
    list["tau_sq"] = tau_sq;
  }
  list["X"] = X;
  if(save_Q){
    list["Q_vals"] = Q_vals;  
  }
  
  if(savepath){
    list["theta_path"] = theta_path;
  }
  return list;
}