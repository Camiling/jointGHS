#' Perform joint GHS with ECM/ICM
#' 
#' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the joint graphical horseshoe
#' 
#' @param X list of the \eqn{K} observed \eqn{n_k} by \eqn{p} data matrices
#' @param theta list of the \eqn{K} initial \eqn{p} by \eqn{p} precision matrices. Optional argument
#' @param sigma list of the \eqn{K} initial \eqn{p} by \eqn{p} covariance matrices. Optional argument
#' @param Lambda_sq list of the \eqn{K} initial \eqn{p} by \eqn{p} matrices of squared local shrinkage parameters. Optional argument
#' @param tau_sq initial values of squared global shrinkage parameters. A vector of length \eqn{K}, or a single value to be used for all networks.
#' @param method the method to use. Default is \code{ECM}. Other options include \code{ICM}
#' @param AIC_selection logical. Should the global shrinkage parameters be selected with AIC? Default \code{TRUE}
#' @param AIC_eps if \code{AIC_selection=TRUE}, the tolerance for the AIC convergence assessment
#' @param tau_sq_min if \code{AIC_selection=TRUE}, the smallest value of \code{tau_sq} to consider  
#' @param tau_sq_stepsize if \code{AIC_selection=TRUE}, the step-size to use in the grid for \code{tau_sq}. Optional argument
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param scale should variables be scaled? Default \code{TRUE}
#' @param verbose logical indicator of printing information at each iteration
#' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for \eqn{p<200}
#' @param save_Q should the value of the objective function at each step be saved?
#' @param fix_tau logical. Should \code{tau_sq} be fixed? Default \code{FALSE}
#' @param boot_check should bootstrapping be performed to check the joint results?
#' @param B if \code{boot_check=TRUE}, the number of bootstrap samples to draw
#' @param nCores if \code{boot_check=TRUE}, how many cores should be used for the bootstrap sampling?
#' 
#' @return a fitted jointGHS object
#' 
#' @importFrom foreach %dopar%
#' 
#' @seealso \code{\link{plot.jointGHS}}, \code{\link{print.jointGHS}} 
#' 
#' @export 
#' 
jointGHS <- function(X, theta=NULL,sigma=NULL,Lambda_sq=NULL, tau_sq = NULL, method= 'ECM', AIC_selection=TRUE, AIC_eps = 1e-1, tau_sq_min =1e-3, tau_sq_stepsize= NULL,
                     epsilon = 1e-5, maxitr = 1e3, scale=TRUE, verbose=TRUE, savepath = FALSE, save_Q = FALSE, fix_tau = FALSE, boot_check=FALSE,
                     B=100,nCores=5){

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
  
  # Select tau for each network using AIC crierion
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
    # Use output to initialize additional parameters
    theta = array(unlist(theta_single), c(p,p,K))
  }
  if(boot_check & AIC_selection){
    perform_boot = c(0,rep(1,B)) # Ensures joint model is run in parallel as well
    doParallel::registerDoParallel(nCores)
    boot.list = foreach (b=1:(B+1)) %dopar% {
      boot_and_joint_parallel_iteration(X, S=S, n.vals = n.vals, p=p, K=K, theta=theta, sigma=sigma, Lambda_sq = Lambda_sq, tau_sq = tau_sq,method=method,
                                        use_ICM = use_ICM, tau_sq_min = tau_sq_min, tau_sq_stepsize = tau_sq_stepsize, epsilon = epsilon, maxitr = maxitr, 
                                        AIC_eps = AIC_eps, perform_boot=perform_boot[b]);
    }
    foreach::registerDoSEQ()
    # Joint results
    out <- boot.list[[1]] 
    # Bootstrap results
    boot.res <- boot.list[-1]
    boot.lambda.list = lapply(boot.res, FUN = function(s) s$Lambda_sq) # A length B list of length K lists
    
  }
  else{
    # Perform joint analysis
    if(AIC_selection) fix_tau = T # If tau has already been selected with the AIC for single networks
    out <- ECM_GHS(S, theta, sigma, Lambda_sq, n.vals, p, K, epsilon, verbose, maxitr, savepath, save_Q, tau_sq, use_ICM = use_ICM, fix_tau = fix_tau)
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
  if(AIC_selection){
    out$Lambda_sq_single = Lambda_sq_single
    out$theta_single = theta_single
    out$sigma_single = sigma_single
    out$AIC_scores = AIC_scores
    out$tau_sq_vals = tau_sq_vals
  }
  if(boot_check){
    out$Lambda_sq_boot =boot.lambda.list
  }

  # Save outputs
  out$tau_sq = tau_sq
  out$epsilon = epsilon
  out$maxitr = maxitr
  class(out) = "jointGHS"
  return(out)
}

#' @keywords internal
boot_and_joint_parallel_iteration = function(X, S, n.vals, p, K, theta, sigma, Lambda_sq, tau_sq, method, use_ICM, tau_sq_min, tau_sq_stepsize, 
                                             epsilon, maxitr, AIC_eps, perform_boot) {
  if(perform_boot==1){
    # Use Bayesian booystrap with weights sampled from the Dirichlet distribution
    res = lapply(1:K, FUN = function(k) fastGHS::fastGHS(X[[k]], theta=theta[,,k], sigma=sigma[,,k], Lambda_sq = Lambda_sq[,,k], tau_sq = tau_sq[k], method = method, 
                           AIC_selection=T, AIC_eps = AIC_eps, tau_sq_min = tau_sq_min, tau_sq_stepsize= tau_sq_stepsize,
                           epsilon = epsilon, maxitr = maxitr, verbose = F, weights=c(gtools::rdirichlet(1, rep(1,nrow(X[[k]])))) ))
    # Save only relevant output
    out = list()
    out$Lambda_sq = lapply(res, FUN = function(s) s$Lambda_sq)
    
  }
  else{
    # Perform the joint analysis
    out = ECM_GHS(S, theta, sigma, Lambda_sq, n.vals, p, K, epsilon, verbose=F, maxitr, savepath=F, save_Q=F, tau_sq, use_ICM = use_ICM, fix_tau = T)
  }
  return(out)
  
}

#' Plot function for S3 class "jointGHS"
#' 
#' This function provides functionality for plotting the results of the jointGHS. Can be used to evaluate the local scales \code{Lambda_sq} of a network from the joint analysis against their Bayesian bootstrap distributions, and to visualize the inferred networks
#' 
#' 
#' @param x An object with S3 class \code{"jointGHS"}
#' @param k the network to plot
#' @param plot_boot should the Bayesian bootstrap sample posterior of some of the Lambda_sq be shown? If \code{FALSE}, the jointGHS networks are visualized instead
#' @param edges the edges for which to show the Bayesian bootstrap posterior distribution. An matrix with \eqn{2} columns. If not specified, \eqn{15} edges are selected at random amongst the infererred jointGHS edge set
#' @param quantiles the quantiles to show when plotting the Bayesian bootstrap posterior distribution of the local scale parameters
#'
#' @seealso \code{\link{jointGHS}}, \code{\link{print.jointGHS}}
#' 
#' @export 
plot.jointGHS = function(x, k=1, plot_boot=TRUE, edges = NA, quantiles = c(0.025, 0.975)){
  gcinfo(FALSE)
  p = ncol(x$theta[[k]])
  if(plot_boot){
    # If no specific edges are requested, sample 16 edges at random from the set of edges identified by jointGHS
    if(any(is.na(edges))){ 
      theta = x$theta[[k]]
      diag(theta) = NA
      edges.all = which(abs(theta)>1e-5,arr.ind=T)
      edges.all = t(apply(edges.all,1,sort))
      edges.all = edges.all[!duplicated(edges.all),]
      if(nrow(edges.all)<16){
        edges = edges.all
      }
      else{
        edges = edges.all[sample(1:nrow(edges.all),15),] 
      }
    }
    if(length(edges)==2){
      edges = t(as.matrix(edges))
    }
    edges = as.matrix(edges)
    if(ncol(edges)!=2 | (is.null(dim(edges)) & length(edges)!=2)){
      stop('Edges must be a matrix with 2 columns, with one row per edge pair. \n')
    }
    if(is.null(x$Lambda_sq_boot)){
      stop('Must provide jointGHS object that was found with option boot_check = TRUE')
    }
    n.edgepairs = nrow(edges) 
    B = length(unlist(lapply(x$Lambda_sq_boot, FUN= function(s) s[[k]][1,2])))
    lambdas = matrix(0,n.edgepairs,B)
    for (e in 1:n.edgepairs){
      lambdas[e,] = unlist(lapply(x$Lambda_sq_boot, FUN= function(s) s[[k]][edges[e,1],edges[e,2]]))
    }
    plot.list = list()
    for(e in 1:n.edgepairs){
      lam.df = data.frame(Lambda_sq=lambdas[e,])
      plot.list[[e]] = ggplot2::ggplot(lam.df,ggplot2::aes(x=Lambda_sq))+ggplot2::labs(title=paste0("Edge (", edges[e,1], ",", edges[e,2],")" ))+
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) + ggplot2::ylab('Frequency')+ggplot2::xlab('Lambda_sq')+
        ggplot2::geom_histogram(ggplot2::aes(y=..density..),color='deepskyblue',fill='deepskyblue',bins=50) + ggplot2::theme(legend.position ="none") + 
        ggplot2::geom_density(alpha=.2, color = 'grey30',fill="deepskyblue")+
        ggplot2::geom_vline(ggplot2::aes(xintercept=x$Lambda_sq[[k]][edges[e,1],edges[e,2]]),color="darkolivegreen3", size=0.5)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=quantile(lambdas[e,],quantiles[1])),color="maroon2", linetype="dashed", size=0.5)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=quantile(lambdas[e,],quantiles[2])),color="maroon2", linetype="dashed", size=0.5)
    }
    # Make grobs object to get legend
    get_legend<-function(quantiles){
      df=data.frame(1:10)
      df[,2] = factor(c(rep('jointGHS estimate',5),rep(paste0(100*abs(diff(quantiles)),'% Credible interval'),5)))
      names(df) = c('x', 'points')
      myggplot = ggplot2::ggplot(df,ggplot2::aes(y=x, group=points)) + ggplot2::geom_line(ggplot2::aes(x=x, color=points, linetype=points)) + 
        ggplot2::scale_linetype_manual(values=c("dashed", "solid"))+ ggplot2::scale_color_manual(values=c('maroon2','darkolivegreen3'))+
        ggplot2::theme(legend.position ="left", legend.title =ggplot2::element_blank())
      tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    legend = get_legend(quantiles)
    plot.list[[(n.edgepairs+1)]] = legend

    dev.new()
    gridExtra::grid.arrange(grobs=plot.list,top= paste0('Bayesian Bootstrap distribution of Lambda_sq for network k=', k))
  }
  else{
    # If bootstrap results are not to be plotted, plot resulting networks
    nets = list()
    for(i in 1:length(x$theta)){
      net = network::network(abs(x$theta[[i]])>1e-5,directed=F)
      nets[[i]] = GGally::ggnet2(net,node.size = 5, edge.size = 0.3,alpha=0.9,mode = "fruchtermanreingold",color = 'dodgerblue')+
                  ggplot2::ggtitle(paste0('Network ',i)) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    }
    dev.new()
    gridExtra::grid.arrange(grobs=nets)
  }
}


#' Print function for S3 class "jointGHS"
#' 
#' This function provides functionality for printing a summary of the results of the jointGHS Bayesian bootstrap check. The local scale parameters from the joint analysis are checked against their Bayesian bootstrap distribution in the single network version
#' 
#' 
#' @param x an object with S3 class \code{"jointGHS"}
#' @param k the network to print the results for
#' @param edges the edges for which to show the Bayesian bootstrap posterior distribution. A matrix with \eqn{2} columns. If not provided, \eqn{16} edges are selected at random amongst the infererred jointGHS edge set
#' @param quantiles the quantiles to show when plotting the Bayesian bootstrap posterior distribution of the local scale parameters
#' @param return_df should the results be returned as a data frame? Default \code{FALSE}
#' 
#' @return object of class \code{"data.frame"}
#'
#' @seealso \code{\link{jointGHS}}, \code{\link{plot.jointGHS}}
#' 
#' @export 
print.jointGHS = function(x, k=1, edges = NA, quantiles = c(0.025, 0.975), return_df=F){
  p = ncol(x$theta[[k]])
  # If no specific edges are requested, sample 16 edges at random from the set of edges identified by jointGHS
  if(any(is.na(edges))){ # Find all edges
    theta = x$theta[[k]]
    diag(theta) = NA
    edges.all = which(abs(theta)>1e-5,arr.ind=T)
    edges.all = t(apply(edges.all,1,sort))
    edges = edges.all[!duplicated(edges.all),]
  }
  if(length(edges)==2){
    edges = t(as.matrix(edges))
  }
  edges = as.matrix(edges)
  if(ncol(edges)!=2 | (is.null(dim(edges)) & length(edges)!=2)){
     stop('Edges must be a matrix with 2 columns, with one row per edge pair. \n')
  }
  if(is.null(x$Lambda_sq_boot)){
    stop('Must provide jointGHS object that was found with option boot_check = TRUE')
  }
  n.edgepairs = nrow(edges) 
  B = length(unlist(lapply(x$Lambda_sq_boot, FUN= function(s) s[[k]][1,2])))
  lambdas = matrix(0,n.edgepairs,B)
  cat(paste0('\nCHECK PARAMETERS OF IDENTIFIED EDGES IN NETWORK k=',k, '\n'))
  cat('==========================================================================================\n')
  cat('Lambda_sq checked against its single-network Bayesian bootstrap distribution \n')
  cat('==========================================================================================\n')
  df = data.frame()
  edges.outside=c()
  for (e in 1:n.edgepairs){
    lambdas[e,] = unlist(lapply(x$Lambda_sq_boot, FUN= function(s) s[[k]][edges[e,1],edges[e,2]]))
    df[e,1] = paste0("(", edges[e,1], ",", edges[e,2],")")
    df[e,2] = round(x$Lambda_sq[[k]][edges[e,1],edges[e,2]],4)
    df[e,3] =  round(quantile(lambdas[e,],quantiles[1]),4)
    df[e,4] =  round(quantile(lambdas[e,],quantiles[2]),4)
    if(df[e,2] < df[e,3] | df[e,2] > df[e,4]){
      df[e,5] = '*'
      edges.outside = c(edges.outside,e)
    }
    else{
      df[e,5] = ''
    }
  }
  names(df) = c('Edge', '\ \ \ jointGHS estimate',paste0('\ \ \ \ \  ',quantiles[1],'-quantile'), paste0('\ \ \ \ \ \ ',quantiles[2],'-quantile'), '')
  print(format(df,justify='left',width=20), right=F)
  cat('==========================================================================================\n')
  cat(paste0('\'*\' Outside ',100*abs(diff(quantiles)),'% credible interval\n\n'))
  cat(paste0('Edge parameters outside intervals: ', length(edges.outside), ' (', 100*round(length(edges.outside)/n.edgepairs,3),'%) \n'))
  cat('Expected number of edge parameters outside intervals: ', round((1-abs(diff(quantiles)))*n.edgepairs), ' (', round(100-100*abs(diff(quantiles)),3),'%) \n\n')

  # Return results as a data frame
  names(df) = c('Edge', 'jointGHS_estimate','quantile1', 'quantile2', 'is_outside')
  df$is_outside = df$is_outside == '*'
  if(return_df==T) return(df)
}




