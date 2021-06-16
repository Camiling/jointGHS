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
#include "E_Nu_GHSlike.h"
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
//' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for p<200
//' @param exist_group logical. Are the variables grouped?
//' @param group grouping information.
//' @param N_groups If exist_group==T, the number of groups
//' @param save_Q should the value of the objective function at each step be saved?
//' @param tau_sq initial value of squared global shrinkage parameter. If exist_group==T, a dummy value should be provided
//' @param Tau_sq if exist_group==T, an \eqn{ngroup} by \eqn{ngroup} matrix of initial values of the squared global shrinkage parameters within and between groups. If exist_group==F, a dummy value should be provided
//' @param machine_eps numerical. The machine precision
//' @param use_ICM logical. Should ICM be used instead of ECM? Default value is false
//' 
//' @return A List with resulting ECM estimates, and saved path and objective function convergence information if requested
//' 
//' @export
// [[Rcpp::export]]
List ECM_GHS(arma::mat X, arma::mat S, arma::mat theta, arma::mat sigma, arma::mat Lambda_sq, double epsilon, bool verbose, int maxitr, bool savepath, int exist_group, arma::uvec group, arma::mat N_groups, bool save_Q, double tau_sq, arma::mat Tau_sq, double machine_eps, bool use_ICM=false, bool fix_tau = false, bool GHS_like = false, bool stop_underflow=false) {

  // get dimensions
  const int M = X.n_cols;
  const int N = X.n_rows;
  
  // For saving variables
  int save_dim;
  if (M < 201 & savepath==true){ // Saving exhausts memory if M>201
    save_dim = maxitr;
  }
  else{
    save_dim = 1;
  }
  arma::cube theta_path(M, M, save_dim);
  arma::vec Q_vals(maxitr);
  double Q_val_old= -numeric_limits<double>::max();
  double Q_val_new;
  
  // initialize intermediate values
  int niter,i, count;
  double eps;
  arma::uvec pseq(M);
  for(i = 0; i < M; i++){
    pseq(i) = i;
  }
  
  // Initialise updates
 arma::mat theta_update = theta; // Save previous estimate to assess convergence
 arma::mat E_Nu_mat(M,M);
 arma::mat E_NuInv(M,M);
 arma::mat E_XiInv(M,M); // Used if variables are grouped
 double E_xiInv; // Used if variables are not grouped
 arma::cube theta_sigma_update(M,M,2);
 eps = epsilon + 1;
 niter = 1;
 count = 0;
 List list;
 
  if(savepath){
    theta_path.slice(0) = theta;
  }
    
  while(eps>epsilon & count < maxitr){
      
    // E-step
    if(GHS_like == false){
      E_NuInv = E_Nu(Lambda_sq);
    }
    else{
      // For the GHS-like prior, we update the matrix Nu of latent shrinkage parameters
      E_Nu_mat = E_Nu_GHSlike(tau_sq, theta);
    }
    
    // The updates below are only relevant for standard GHS
    if (exist_group>0){
      E_XiInv = E_xi_group(Tau_sq);
    }
    else if (GHS_like == false){
      E_xiInv = E_xi(tau_sq);
    }
    // Trick for simplicity of computations of ICM: multiply expectations by 2 to get modes
    if(use_ICM){
      if(exist_group>0){
        E_XiInv = 2*E_XiInv;
      }
      else{
        E_xiInv = 2*E_xiInv; 
      }
      E_NuInv = 2*E_NuInv;
    }
    
    // M-step
    if (exist_group>0){
      Lambda_sq = M_lambda(N, M, theta, E_NuInv, exist_group, group, Tau_sq, machine_eps, stop_underflow);
      if (fix_tau == false){
        Tau_sq = M_tau_group(M, theta, Lambda_sq, exist_group, group, N_groups, E_XiInv, machine_eps, stop_underflow); 
      }
      theta_sigma_update = M_theta(N, M, theta, S, sigma, Lambda_sq, pseq, exist_group, group, Tau_sq, machine_eps, stop_underflow);
    }
    else{
      // For standard GHS, we update Lambda_sq, tau and theta in the M-step
      if (GHS_like == false){
        Lambda_sq = M_lambda(N, M, theta, E_NuInv, exist_group, group, S, machine_eps, stop_underflow, tau_sq);
        if (fix_tau == false){
          tau_sq = M_tau(M, theta, Lambda_sq, E_xiInv, machine_eps, stop_underflow);
        }
        theta_sigma_update = M_theta(N, M, theta, S, sigma, Lambda_sq, pseq, exist_group, group, S,machine_eps, stop_underflow, tau_sq); // Pass S as dummy argument
      }
      // For the GHS-like prior, we update theta only in the M-step
      else {
        // Use trick to reuse code: 
        theta_sigma_update = M_theta(N, M, theta, S, sigma, E_Nu_mat, pseq, exist_group, group, S,machine_eps, stop_underflow, tau_sq/2); // Pass S as dummy argument
      }
    }

    theta_update = theta_sigma_update.slice(0);
    sigma = theta_sigma_update.slice(1);

    if(savepath){
      theta_path.slice(count+1) = theta_update;
    }
    // Evaluate objective function if ICM is used, or if it is to be saved
    if((use_ICM || save_Q) & (exist_group>0)){
      Q_val_new = Q_val(N, M, theta, S, Lambda_sq, E_NuInv, exist_group, group, Tau_sq, E_XiInv); 
    }else if (use_ICM || save_Q) {
      Q_val_new = Q_val(N, M, theta, S, Lambda_sq, E_NuInv, exist_group, group, S, S, tau_sq, E_xiInv);  // Pass S as dummy argument
    }
    // If ICM is used, check that objective function value has increased
    if(use_ICM & (Q_val_old>Q_val_new)){
      Rcout << "Error: objective function decreasing" << endl;
      list["Q_val_old"] = Q_val_old;
      list["Q_val_new"] = Q_val_new;
      return list;
    }
    if(save_Q){
      Q_vals(count) = Q_val_new;
    }
    // Update estimate
    if (use_ICM || save_Q) {
      Q_val_old = Q_val_new;
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
  if(save_Q){
    Q_vals= Q_vals.head(count);
  }
  // Save results
  list["S"] = S;
  list["theta"] = theta;  
  list["sigma"] = sigma;
  list["X"] = X;
  if(GHS_like == false){
    list["Lambda_sq"] = Lambda_sq;
    if(exist_group>0){
      list["tau_sq"] = Tau_sq;
    }
    else{
      list["tau_sq"] = tau_sq;
    }
    if(use_ICM){
      list["Nu"] = 1/E_NuInv;
      if(exist_group>0){
        list["Xi"] = 1/E_XiInv;
      }
      else{
        list["xi"] = 1/E_xiInv;
      }
    }
  } 
  else{
    list["a"] = tau_sq;
    list["Nu"] = E_Nu_mat;
  }
  if(save_Q){
    list["Q_vals"] = Q_vals;  
  }
  
  if(savepath){
    list["theta_path"] = theta_path;
  }
  return list;
}