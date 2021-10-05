#' Perform joint GHS with ECM/ICM
#' 
#' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the joint graphical horseshoe
#' 
#' @param X list of the \eqn{K} observed \eqn{n_k} by \eqn{p} data matrices
#' @param theta list of the \eqn{K} initial \eqn{p} by \eqn{p} precision matrices. Optional argument
#' @param sigma list of the \eqn{K} initial \eqn{p} by \eqn{p} covariance matrices. Optional argument
#' @param Lambda_sq list of the \eqn{K} initial \eqn{p} by \eqn{p} matrices of squared local shrinkage parameters.
#' @param tau_sq initial values of squared global shrinkage parameters. A vector of length \eqn{K}
#' @param method the method to use. Default is \eqn{ECM}. Other options include \eqn{ICM}
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param verbose logical indicator of printing information at each iteration
#' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for p<200
#' @param save_Q should the value of the objective function at each step be saved?
#' @param fix_tau logical. Should tau be fixed?
#' @return a fitted EMGS object
#' @export 
#' 
jointGHS <- function(X, theta=NULL,sigma=NULL,Lambda_sq=NULL, tau_sq = NULL, method= 'ECM',epsilon = 1e-5, maxitr = 1e5, verbose=TRUE, savepath = FALSE, save_Q = FALSE, fix_tau = TRUE){

  p <- dim(X[[1]])[2]
  K <- length(X)
  
  for(k in 1:K){
    X[[k]] = as.matrix(X[[k]])
  }
  # Assign initial values, unless provided
  # Create arrays
  if(is.null(theta)){
    theta <- array(0,dim=c(p,p,K))
    for(k in 1:K){
      theta[,,k] <- diag(1,p)
    }
  }
  else{
    theta_list = theta
    theta = array(0,dim=c(p,p,K))
    for(k in 1:K){
      theta[,,k] = as.matrix(theta[[k]])
      if(!isSymmetric(theta[,,k])){
        theta[,,k] = as.matrix(Matrix::forceSymmetric(theta[,,k]))
        cat('Initial theta of data set ', k , ' not symmetric, forcing symmetry...\n')
      }
      if(ncol(theta[,,k])!= p | nrow(theta[,,k])!=p | !matrixcalc::is.positive.definite(theta[,,k]+0)){
        cat('Error: initial theta of data set ', k, ' must be pxp and positive definite \n')
        return()
      } 
    }

  }
  if(is.null(sigma)){
    sigma = array(0,dim=c(p,p,K))
    for(k in 1:K){
      sigma[,,k] <- diag(1,p)
    }
  }
  else{
    sigma_list = sigma
    sigma = array(0,dim=c(p,p,K))
    for(k in 1:K){
      sigma[,,k] = as.matrix(sigma_list[[k]])
      if(!isSymmetric(sigma[,,k])){
        sigma[,,k]=as.matrix(Matrix::forceSymmetric(sigma[,,k]))
        cat('Initial sigma of data set ', k, ' not symmetric, forcing symmetry...\n')
      }
      if(ncol(sigma[,,k])!= p | nrow(sigma[,,k] )!=p | !matrixcalc::is.positive.definite(as.matrix(sigma[,,k]+0))){
        cat('Error: initial sigma of data set ', k, ' must be pxp and positive definite \n')
        return()
      } 
    }
  }
  if(is.null(Lambda_sq)){
    Lambda_sq = array(1,dim=c(p,p,K))
  }
  else{
    Lambda_sq_list = Lambda_sq
    Lambda_sq = array(0,dim=c(p,p,K))
    for(k in 1:K){
      Lambda_sq[,,k] = Lambda_sq_list[[k]]
      if(any(Lambda_sq[,,k]<0)){
        cat('Error: negative Lambda_sq values are not allowed \n')
        return()
      } 
    }

  }

  n.vals <- unlist(lapply(X, nrow))
  S_list <- lapply(X, FUN = function(x) t(x) %*% x)
  S = array(0,dim=c(p,p,K))
  for(k in 1:K){
    S[,,k] = S_list[[k]] 
  }
  
  if(is.null(tau_sq)){
    tau_sq <- rep(1, K)
  }
  else{
    if(length(tau_sq)!=K){
      cat('Error: number of tau_sq values must match number of networks \n')
      return()
    }
       
    if(any(tau_sq<0)){
      cat('Error: negative tau_sq is not allowed \n')
      return()
    }
  }

  if(method=='ECM'){
    out <- ECM_GHS(S, theta , sigma, Lambda_sq, n.vals, p, K, epsilon, verbose, maxitr, savepath, save_Q, tau_sq, use_ICM = FALSE, fix_tau = fix_tau)
  }
  else if(method=='ICM'){
    out <- ECM_GHS(S, theta , sigma, Lambda_sq, n.vals, p, K, epsilon, verbose, maxitr, savepath, save_Q, tau_sq, use_ICM = TRUE, fix_tau = fix_tau)
  }
  else {
    out <- NULL
  }
  
  # Make lists for the output
  
  S_list =  out$S
  theta_list =  out$theta
  sigma_list =  out$sigma
  Lambda_sq_list = out$Lambda_sq
  out$S = list()
  out$theta = list()
  out$sigma = list()
  out$Lambda_sq = list()
  for (k in 1:K){
    out$S[[k]] = S_list[,,k]
    out$theta[[k]] = theta_list[,,k]
    out$sigma[[k]] = sigma_list[,,k]
    out$Lambda_sq[[k]] = Lambda_sq_list[,,k]
  }

  # Save outputs
  out$epsilon = epsilon
  out$maxitr = maxitr
  class(out) = "fastGHS"
  return(out)
}
