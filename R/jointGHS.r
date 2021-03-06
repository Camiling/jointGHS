#' Perform joint GHS with ECM/ICM
#' 
#' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the joint graphical horseshoe
#' 
#' @param X list of the \eqn{K} observed \eqn{n_k} by \eqn{p} data matrices
#' @param epsilon tolerance for the convergence assessment
#' @param AIC_eps if \code{AIC_selection=TRUE}, the tolerance for the AIC convergence assessment
#' @param theta list of the \eqn{K} initial \eqn{p} by \eqn{p} precision matrices. Optional argument
#' @param sigma list of the \eqn{K} initial \eqn{p} by \eqn{p} covariance matrices. Optional argument
#' @param Lambda_sq list of the \eqn{K} initial \eqn{p} by \eqn{p} matrices of squared local shrinkage parameters. Optional argument
#' @param tau_sq initial values of squared global shrinkage parameters. A vector of length \eqn{K}, or a single value to be used for all networks.
#' @param tau_sq_min if \code{AIC_selection=TRUE}, the smallest value of \code{tau_sq} to consider  
#' @param tau_sq_stepsize if \code{AIC_selection=TRUE}, the step-size to use in the grid for \code{tau_sq}. Optional argument
#' @param tau_sq_max if \code{AIC_selection=TRUE}, the largest value of \code{tau_sq} to consider 
#' @param stop_overflow should measures be taken to avoid overflow? Only necessary for very dense graphs with many large nonzero elements. 
#' @param maxitr maximum number of iterations
#' @param scale should variables be scaled? Default \code{TRUE}
#' @param verbose logical indicator of printing information at each iteration
#' @param fix_tau logical. Should \code{tau_sq} be fixed? Default \code{FALSE}
#' @param AIC_selection logical. Should the global shrinkage parameters be selected with AIC? Default \code{TRUE}
#' @param boot_check should bootstrapping be performed to check the joint results?
#' @param B if \code{boot_check=TRUE}, the number of bootstrap samples to draw
#' @param nCores if \code{boot_check=TRUE}, how many cores should be used for the bootstrap sampling?
#' 
#' @return a fitted \code{jointGHS} object
#' 
#' @importFrom foreach %dopar%
#' @importFrom ggplot2 ..density..
#' 
#' @seealso \code{\link{plot_jointGHS}}, \code{\link{print_jointGHS}} 
#' 
#' @export 
#' 
jointGHS <- function(X, epsilon = 1e-5, AIC_eps = 1e-1, theta=NULL,sigma=NULL,Lambda_sq=NULL, tau_sq = NULL, tau_sq_min =1e-3, tau_sq_stepsize= NULL,tau_sq_max=20,
                     stop_overflow=FALSE, maxitr = 1e3, scale=TRUE, verbose=TRUE, fix_tau = FALSE, AIC_selection=TRUE, boot_check=FALSE, B=100,nCores=5){

  # If only one data set is provided, single-network version is used instead.
  if(!is.list(X)){
    if(is.list(theta) | is.list(sigma) | is.list(Lambda_sq) | length(tau_sq)!=1){
      stop('Cannot provide several starting values for single network analysis.')
    }
    # All other checks are performed within fastGHS
    cat('Single data set provided: performing single-network ECM for GHS \n')
    out <- fastGHS::fastGHS(X, theta=theta, sigma=sigma, Lambda_sq = Lambda_sq, tau_sq = tau_sq, AIC_selection = AIC_selection, AIC_eps = AIC_eps,
                            tau_sq_min=tau_sq_min, tau_sq_stepsize = tau_sq_stepsize, tau_sq_max=tau_sq_max,
                            epsilon = epsilon, maxitr = maxitr, scale=scale, verbose = verbose, fix_tau = fix_tau, stop_overflow = stop_overflow)
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
    if(any(is.null(tau_sq))){
      stop('Cannot select fix_tau=T without providing a tau_sq value \n')
    }
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

  # Select tau for each network using AIC crierion
  if(AIC_selection & !fix_tau){
    doParallel::registerDoParallel(K)
    res.list = foreach::foreach (k=1:K) %dopar% {
      fastGHS::fastGHS(X[[k]], theta=theta[,,k], sigma=sigma[,,k], Lambda_sq = Lambda_sq[,,k], tau_sq = tau_sq[k],
                       AIC_selection=AIC_selection, AIC_eps = AIC_eps, tau_sq_min =tau_sq_min, tau_sq_stepsize= tau_sq_stepsize,tau_sq_max=tau_sq_max,
                       epsilon = epsilon, maxitr = maxitr, verbose = F, stop_overflow = stop_overflow);
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
  if(boot_check & AIC_selection & !fix_tau){
    perform_boot = c(0,rep(1,B)) # Ensures joint model is run in parallel as well
    doParallel::registerDoParallel(nCores)
    boot.list = foreach::foreach (b=1:(B+1)) %dopar% {
      boot_and_joint_parallel_iteration(X, S=S, n.vals = n.vals, p=p, K=K, theta=theta, sigma=sigma, Lambda_sq = Lambda_sq, tau_sq = tau_sq,
                                        tau_sq_min = tau_sq_min, tau_sq_stepsize = tau_sq_stepsize, tau_sq_max = tau_sq_max, epsilon = epsilon, maxitr = maxitr, 
                                        AIC_eps = AIC_eps, perform_boot=perform_boot[b],stop_overflow=stop_overflow);
    }
    foreach::registerDoSEQ()
    # Joint results
    out <- boot.list[[1]] 
    # Bootstrap results
    boot.res <- boot.list[-1]
    boot.theta.list = lapply(boot.res, FUN = function(s) s$theta) # A length B list of length K lists
    
  }
  else{
    # Perform joint analysis
    if(AIC_selection) fix_tau = T # If tau has already been selected with the AIC for single networks
    out <- ECM_GHS(S, theta, sigma, Lambda_sq, n.vals, p, K, epsilon, verbose, maxitr, tau_sq, fix_tau = fix_tau)
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
  if(AIC_selection & exists('AIC_scores')){
    out$Lambda_sq_single = Lambda_sq_single
    out$theta_single = theta_single
    out$sigma_single = sigma_single
    out$AIC_scores = AIC_scores
    out$tau_sq_vals = tau_sq_vals
  }
  if(boot_check & exists('boot.theta.list')){
    out$theta_boot =boot.theta.list
  }

  # Save outputs
  out$tau_sq = tau_sq
  out$epsilon = epsilon
  out$maxitr = maxitr
  class(out) = "jointGHS"
  return(out)
}

#' Perform bootstrap iterarion
#' 
#' This function performs one Bayesian bootstrap iteration
#' 
#' @param X list of the \eqn{K} observed \eqn{n_k} by \eqn{p} data matrices
#' @param S list of the \eqn{p} by \eqn{p} scatter matrices of the \eqn{K} data sets
#' @param n.vals vector of length \eqn{K} of the number of obervations in each data set.
#' @param p number of nodes
#' @param K number of networks
#' @param theta list of the \eqn{K} initial \eqn{p} by \eqn{p} precision matrices. Optional argument
#' @param sigma list of the \eqn{K} initial \eqn{p} by \eqn{p} covariance matrices. Optional argument
#' @param Lambda_sq list of the \eqn{K} initial \eqn{p} by \eqn{p} matrices of squared local shrinkage parameters. Optional argument
#' @param tau_sq initial values of squared global shrinkage parameters. A vector of length \eqn{K}, or a single value to be used for all networks.
#' @param tau_sq_min if \code{AIC_selection=TRUE}, the smallest value of \code{tau_sq} to consider  
#' @param tau_sq_stepsize if \code{AIC_selection=TRUE}, the step-size to use in the grid for \code{tau_sq}. Optional argument
#' @param tau_sq_max if \code{AIC_selection=TRUE}, the largest value of \code{tau_sq} to consider  
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param AIC_eps if \code{AIC_selection=TRUE}, the tolerance for the AIC convergence assessment
#' @param perform_boot should bootstrapping be performed? if 1, yes, else perform joint
#' @param stop_overflow should measures be taken to avoid overflow? Only necessary for very dense graphs with many large nonzero elements.
#' 
#' 
#' @keywords internal
boot_and_joint_parallel_iteration = function(X, S, n.vals, p, K, theta, sigma, Lambda_sq, tau_sq, tau_sq_min, tau_sq_stepsize, tau_sq_max,
                                             epsilon, maxitr, AIC_eps, perform_boot,stop_overflow) {
  if(perform_boot==1){
    # Use Bayesian bootstrap with weights sampled from the Dirichlet distribution
    res = lapply(1:K, FUN = function(k) fastGHS::fastGHS(X[[k]], theta=theta[,,k], sigma=sigma[,,k], Lambda_sq = Lambda_sq[,,k], tau_sq = tau_sq[k], 
                           AIC_selection=T, AIC_eps = AIC_eps, tau_sq_min = tau_sq_min, tau_sq_stepsize= tau_sq_stepsize, tau_sq_max = tau_sq_max,
                           epsilon = epsilon, maxitr = maxitr, verbose = F, weights=c(gtools::rdirichlet(1, rep(1,nrow(X[[k]])))), stop_overflow = stop_overflow ))
    # Save only relevant output
    out = list()
    out$theta = lapply(res, FUN = function(s) s$theta)
    
  }
  else{
    # Perform the joint analysis
    out = ECM_GHS(S, theta, sigma, Lambda_sq, n.vals, p, K, epsilon, verbose=F, maxitr, tau_sq, fix_tau = T)
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
#' @param percentile the percentile to show when plotting the Bayesian bootstrap posterior distribution of the local scale parameters.
#' @param quantiles if provided, the quantiles to show when plotting the Bayesian bootstrap posterior distribution of the local scale parameters. Overrides \code{percentile} argument
#' @param show_single should the estimate from the single network analysis be included as a bar?
#' @param true_theta if provided, the 'true' value of \code{theta} to include as a bar in the plot
#'
#' @seealso \code{\link{jointGHS}}, \code{\link{print_jointGHS}}
#' 
#' @importFrom ggplot2 ..density..
#' 
#' @export 
plot_jointGHS = function(x, k=1, plot_boot=TRUE, edges = NA, percentile=95, quantiles = NULL, show_single=FALSE, true_theta=NULL){
  p = ncol(x$theta[[k]])
  quantiles=c(0,percentile/100) # Show 0-quantile as well as percentile/100-quantile
  if(plot_boot){
    # If no specific edges are requested, sample 15 edges at random from the set of edges identified by jointGHS
    if(any(is.na(edges))){ 
      theta = x$theta[[k]]
      diag(theta) = NA
      edges.all = which(abs(theta)>1e-5,arr.ind=T)
      edges.all = t(apply(edges.all,1,sort))
      edges.all = edges.all[!duplicated(edges.all),]
      if(any(is.na(edges.all)) | any(is.null(edges.all)) | length(edges.all) ==0){
        stop('No edges in network to print')
      }
      if(length(edges.all) == 2){
        edges = edges.all
      }
      else if(nrow(edges.all)<15){
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
    
    # If we are checking theta
    x$Lambda_sq_boot = x$theta_boot
    x$Lambda_sq = lapply(x$theta, cov2cor) # Get correlation so we can compare regardless of diagonal element values
    x$Lambda_sq_single = lapply(x$theta_single, cov2cor)
    if(!is.null(true_theta)){
      true_theta = cov2cor(true_theta)
    }
    
    if(is.null(x$Lambda_sq_boot)){
      stop('Must provide jointGHS object that was found with option boot_check = TRUE')
    }
    n.edgepairs = nrow(edges) 
    B = length(unlist(lapply(x$Lambda_sq_boot, FUN= function(s) s[[k]][1,2])))
    lambdas = matrix(0,n.edgepairs,B)
    for (e in 1:n.edgepairs){
      lambdas[e,] = unlist(lapply(x$Lambda_sq_boot, FUN= function(s) cov2cor(s[[k]])[edges[e,1],edges[e,2]] ))
    }
    xlab.name = 'theta'
    plot.list= lapply(1:n.edgepairs, FUN = function(j) ggplot2::ggplot(data.frame(Lambda_sq=lambdas[j,]),ggplot2::aes(x=Lambda_sq))+ggplot2::labs(title=paste0("Edge (", edges[j,1], ",", edges[j,2],")" ))+
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) + ggplot2::ylab('Frequency')+ggplot2::xlab(xlab.name)+
        ggplot2::geom_histogram(ggplot2::aes(y=..density..),color='deepskyblue',fill='deepskyblue',bins=50) + ggplot2::theme(legend.position ="none") + 
        ggplot2::geom_density(alpha=.2, color = 'grey30',fill="deepskyblue")+
        ggplot2::geom_vline(ggplot2::aes(xintercept=x$Lambda_sq[[k]][edges[j,1],edges[j,2]]),color="darkolivegreen3", size=0.5)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=sign(mean(sign(lambdas[j,])))*quantile(abs(lambdas[j,]),quantiles[2])),color="maroon2", linetype="dashed", size=0.5)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=sign(mean(sign(lambdas[j,])))*quantile(abs(lambdas[j,]),quantiles[1])),color="maroon2", linetype="dashed", size=0.5)+
        {if(show_single) ggplot2::geom_vline(ggplot2::aes(xintercept=x$Lambda_sq_single[[k]][edges[j,1],edges[j,2]]),color="darkorange", linetype="solid", size=0.5)}+
        {if(!is.null(true_theta)) ggplot2::geom_vline(ggplot2::aes(xintercept=true_theta[edges[j,1],edges[j,2]]),color="aquamarine", linetype="solid", size=0.5)})
    
    # Make grobs object to get legend
    get_legend<-function(quantiles, show_single=F, true_theta=NULL){
      df=data.frame(1:10)
      if(show_single & !is.null(true_theta)){
        df = data.frame(1:20)
      }
      else if(show_single | !is.null(true_theta)){
        df = data.frame(1:15)
      }
      names.all = c(rep('jointGHS estimate',5),rep(paste0(quantiles[1], ' and ', quantiles[2]*100 ,'% empirical percentiles'),5))
      if(show_single) names.all = c(names.all, rep('single estimate',5))
      if(!is.null(true_theta)) names.all =c(names.all, rep('truth',5))
      df[,2]= factor(names.all)
      
      names(df) = c('x', 'points')
      col.vals = c('maroon2','darkolivegreen3')
      if(show_single) col.vals = c(col.vals,'darkorange')
      if(!is.null(true_theta)) col.vals = c(col.vals,'aquamarine')
      type.vals = c("dashed", "solid")
      if(show_single) type.vals = c(type.vals,'solid')
      if(!is.null(true_theta)) type.vals = c(type.vals,'solid')
      
      myggplot = ggplot2::ggplot(df,ggplot2::aes(y=x, group=points)) + ggplot2::geom_line(ggplot2::aes(x=x, color=points, linetype=points)) + 
        ggplot2::scale_linetype_manual(values=type.vals)+ ggplot2::scale_color_manual(values=col.vals)+
        ggplot2::theme(legend.position ="left", legend.title =ggplot2::element_blank())

      tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    legend = get_legend(quantiles, show_single, true_theta)
    plot.list[[(n.edgepairs+1)]] = legend

    dev.new()
    gridExtra::grid.arrange(grobs=plot.list,top= paste0('Bayesian Bootstrap distribution of ', xlab.name, ' for network k=', k))
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
#' @param edges the edges for which to show the Bayesian bootstrap posterior distribution. A matrix with \eqn{2} columns, or alternatively the string \code{'all'}. If not provided, \eqn{16} edges are selected at random amongst the infererred jointGHS edge set
#' @param percentile the empirical percentile of the Bayesian bootstrap posterior distribution of the local scale parameters we want to show. Can also be a vector of length \eqn{2} for two-sided intervals
#' @param return_df should the results be returned as a data frame? Default \code{FALSE}
#'
#' @return object of class \code{"data.frame"}
#'
#' @seealso \code{\link{jointGHS}}, \code{\link{plot_jointGHS}}
#' 
#' @export 
print_jointGHS = function(x, k=1, edges = NA, percentile = 95, return_df=F){
  p = ncol(x$theta[[k]])
  if(length(percentile)==2){
    quantiles=c(percentile[1]/100, percentile[2]/100)
    onesided=F
  }
  else if(length(percentile)==1){
    quantiles = percentile/100
    onesided=T
  }
  else{
    stop('percentile must be of length 1 or 2')
  }
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
  if(any(edges=='all')){
    edges = t(combn(1:p, 2))
  }
  edges = as.matrix(edges)
  if(ncol(edges)!=2 | (is.null(dim(edges)) & length(edges)!=2)){
     stop('Edges must be a matrix with 2 columns, with one row per edge pair. \n')
  }
  # Checking theta
  x$Lambda_sq_boot = x$theta_boot
  x$Lambda_sq = x$theta
  if(is.null(x$Lambda_sq_boot)){
    stop('Must provide jointGHS object that was found with option boot_check = TRUE')
  }
  n.edgepairs = nrow(edges) 
  B = length(unlist(lapply(x$Lambda_sq_boot, FUN= function(s) s[[k]][1,2])))
  lambdas = matrix(0,n.edgepairs,B)
  var.name = 'theta'
  cat(paste0('\nCHECK PARAMETERS OF IDENTIFIED EDGES IN NETWORK k=',k, '\n'))
  if(onesided) partbar = '==========================================================================================\n'
  else partbar = '======================================================================================================================\n'
  cat(partbar)
  cat(paste0(var.name,' checked against its single-network Bayesian bootstrap distribution \n'))
  cat(partbar)
  df = data.frame()
  edges.outside=c()
  for (e in 1:n.edgepairs){
    lambdas[e,] = unlist(lapply(x$Lambda_sq_boot, FUN= function(s) cov2cor(s[[k]])[edges[e,1],edges[e,2]] ))
    joint.est =  cov2cor(x$Lambda_sq[[k]])[edges[e,1],edges[e,2]]
    df[e,1] = paste0("(", edges[e,1], ",", edges[e,2],")")
    df[e,2] = round(joint.est,4)
    if(onesided){
      quant1 =  sign(mean(sign(lambdas[e,])))*quantile(abs(lambdas[e,]),quantiles) # Get the right sign
      df[e,3] =  round(quant1,4)
      if(abs(df[e,2]) > abs(df[e,3]) | sign(df[e,2]) != sign(df[e,3])){ # Must have the same sign
        df[e,4] = '*'
        edges.outside = c(edges.outside,e)
      }
      else{
        df[e,4] = ''
      }
    }
    else{
      quant1 = sign(mean(sign(lambdas[e,])))*quantile(abs(lambdas[e,]),quantiles[1])
      quant2 = sign(mean(sign(lambdas[e,])))*quantile(abs(lambdas[e,]),quantiles[2])
      df[e,3] =  round(quant1,4)
      df[e,4] =  round(quant2,4)
      if(abs(df[e,2]) < abs(df[e,3]) | abs(df[e,2]) > abs(df[e,4]) | !(sign(df[e,2]) %in% c(sign(df[e,3]), sign(df[e,4]))) ){ # Must have the same sign too
        df[e,5] = '*'
        edges.outside = c(edges.outside,e)
      }
      else{
        df[e,5] = ''
      }
    }

  }
  if(onesided){
    names(df) = c('Edge', '\ \ \ \ \ \ \ \ \  jointGHS estimate',paste0('\ \ \ ',quantiles*100,'% empirical percentile'),  '')
  }
  else{
    names(df) = c('Edge', '\ \ \ \ \ \ \ \ \  jointGHS estimate',paste0('\ \ \  ',quantiles[1]*100,'% empirical percentile'), paste0('\ \ \ ',quantiles[2]*100,'% empirical percentile'), '')
  }
  print(format(df,justify='left',width=27), right=F)
  cat(partbar)
  if(onesided) { 
    cat(paste0('\'*\' Larger than ',100*quantiles,'% empirical percentile in absolute value \n\n')) 
    cat(paste0('Edge parameters larger than ',100*quantiles,'% empirical percentile in absolute value: ', length(edges.outside), ' (', 100*round(length(edges.outside)/n.edgepairs,3),'%) \n'))
    cat('Expected number of edge parameters larger than ', 100*quantiles,'% empirical percentile in absolute value:', round((1-quantiles)*n.edgepairs), ' (', round(100-100*quantiles,3),'%) \n\n')
    names(df) = c('Edge', 'jointGHS_estimate','percentile', 'is_outside')
  }
  else {
    cat(paste0('\'*\' Outside ',100*abs(diff(quantiles)),'% empirical interval \n\n'))
    cat(paste0('Edge parameters outside interval: ', length(edges.outside), ' (', 100*round(length(edges.outside)/n.edgepairs,3),'%) \n'))
    cat('Expected number of edge parameters outside interval: ', round((1-abs(diff(quantiles)))*n.edgepairs), ' (', round(100-100*abs(diff(quantiles)),3),'%) \n\n')
    names(df) = c('Edge', 'jointGHS_estimate','percentile1', 'percentile2', 'is_outside')
  }
  # Return results as a data frame
  df$is_outside = df$is_outside == '*'
  if(return_df==T) return(df)
}




