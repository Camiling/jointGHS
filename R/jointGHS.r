#' Perform joint GHS with ECM/ICM
#' 
#' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the joint graphical horseshoe
#' 
#' @param X list of the \eqn{K} observed \eqn{n_k} by \eqn{p} data matrices
#' @param theta list of the \eqn{K} initial \eqn{p} by \eqn{p} precision matrices. Optional argument
#' @param sigma list of the \eqn{K} initial \eqn{p} by \eqn{p} covariance matrices. Optional argument
#' @param Lambda_sq list of the \eqn{K} initial \eqn{p} by \eqn{p} matrices of squared local shrinkage parameters. Optional argument
#' @param tau_sq initial values of squared global shrinkage parameters. A vector of length \eqn{K}, or a single value to be used for all networks
#' @param method the method to use. Default is \eqn{ECM}. Other options include \eqn{ICM}
#' @param AIC_selection logical. Should the global shrinkage parameters be selected with AIC?
#' @param AIC_eps if AIC_selection == TRUE, the convergence tolerance for the AIC convergence assessment
#' @param tau_sq_min if AIC_selection == TRUE, the smallest value of tau_sq to consider  
#' @param tau_sq_stepsize if AIC_selection == TRUE, the step-size to use in the grid for tau_sq. Optional argument
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param scale should variables be scaled?
#' @param verbose logical indicator of printing information at each iteration
#' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for p<200
#' @param save_Q should the value of the objective function at each step be saved?
#' @param fix_tau logical. Should tau be fixed?
#' @return a fitted EMGS object
#' 
#' @importFrom foreach %dopar%
#' 
#' @export 
#' 
jointGHS <- function(X, theta=NULL,sigma=NULL,Lambda_sq=NULL, tau_sq = NULL, method= 'ECM', AIC_selection=TRUE, AIC_eps = 1e-1, tau_sq_min =1e-3, tau_sq_stepsize= NULL,
                     epsilon = 1e-5, maxitr = 1e3, scale=FALSE, verbose=TRUE, savepath = FALSE, save_Q = FALSE, fix_tau = FALSE){

  # If only one data set is provided, single-network version is used instead.
  if(!is.list(X)){
    if(is.list(theta) | is.list(sigma) | is.list(Lambda_sq) | length(tau_sq)!=1){
      stop('Cannot provide several starting values for single network analysis.')
    }
    # All other checks are performed within fastGHS
    cat('Single data set provided: performing single-network ECM for GHS \n')
    out <- fastGHS::fastGHS(X, theta=theta, sigma=sigma, Lambda_sq = Lambda_sq, tau_sq = tau_sq, method = method, 
                            epsilon = epsilon, maxitr = maxitr, verbose = verbose, savepath=savepath,save_Q=save_Q, fix_tau = fix_tau)
    class(out) = "jointGHS"
    return(out)
  }
  
  p <- dim(X[[1]])[2]
  K <- length(X)
  
  for(k in 1:K){
    X[[k]] = as.matrix(X[[k]])
    if(scale){
      X[[k]] = scale(X[[k]])
    }
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
        stop('Initial theta of data set ', k, ' must be pxp and positive definite \n')
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
        stop('Initial sigma of data set ', k, ' must be pxp and positive definite \n')
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
        stop('Negative Lambda_sq values are not allowed \n')
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
    if(length(tau_sq)==1){
      tau_sq = rep(tau_sq, K)
    }
    else if(length(tau_sq)!=K){
      stop('Number of tau_sq values must match number of networks, or a single value must be provided for all. \n')
    }
       
    if(any(tau_sq<0)){
      stop('Negative tau_sq is not allowed \n')
    }
  }
  
  if(!method %in% c('ECM', 'ICM') ){
    stop('Method must be either ECM or ICM')
  }
  use_ICM = method=='ICM'
  
  # Select tau for each network using AIC crieria
  if(AIC_selection){
    doParallel::registerDoParallel(K)
    res.list = foreach (k=1:K) %dopar% {
      fastGHS::fastGHS(X[[k]], theta=theta[,,k], sigma=sigma[,,k], Lambda_sq = Lambda_sq[,,k], tau_sq = tau_sq[k], method = method, 
                       AIC_selection=AIC_selection, AIC_eps = AIC_eps, tau_sq_min =tau_sq_min, tau_sq_stepsize= tau_sq_stepsize,
                       epsilon = epsilon, maxitr = maxitr, verbose = F);
    }
    foreach::registerDoSEQ()
    # New tau_sq values
    tau_sq = unlist(lapply(res.list, FUN = function(s) s$tau_sq))
    # List of vectors of AIC score trajectories
    AIC_scores = lapply(res.list, FUN = function(s) s$AIC_scores)
    # Save single-network output
    Lambda_sq_single = lapply(res.list, FUN = function(s) s$Lambda_sq)
    theta_single = lapply(res.list, FUN = function(s) s$theta)
    sigma_single = lapply(res.list, FUN = function(s) s$sigma)
    tau_sq_vals = lapply(res.list, FUN = function(s) s$tau_sq_vals)
    # Use output to initialise 
    theta = array(unlist(theta_single), c(p,p,K))
    Lambda_sq = array(unlist(Lambda_sq_single), c(p,p,K))
    sigma = array(unlist(sigma_single), c(p,p,K))
    
  }
  
  # Perform joint analysis
  if(AIC_selection) fix_tau = T
  out <- ECM_GHS(S, theta, sigma, Lambda_sq, n.vals, p, K, epsilon, verbose, maxitr, savepath, save_Q, tau_sq, use_ICM = use_ICM, fix_tau = fix_tau)

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
  if(AIC_selection){
    out$Lambda_sq_single = Lambda_sq_single
    out$theta_single = theta_single
    out$sigma_single = sigma_single
    out$AIC_scores = AIC_scores
    out$tau_sq_vals = tau_sq_vals
  }

  # Save outputs
  out$tau_sq = tau_sq
  out$epsilon = epsilon
  out$maxitr = maxitr
  class(out) = "jointGHS"
  return(out)
}
