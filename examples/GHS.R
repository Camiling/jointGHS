# Standard MCMC for graphical horseshoe

GHS =  function(S,n,burnin=500,nmc=3000){
  # GHS MCMC sampler using data-augmented
  # block (column-wise) Gibbs sampler
  # Input:
  #     S = Y'*Y : sample covariance matrix * n
  #     n: sample size
  #     burnin, nmc : number of MCMC burnins and saved samples
  
  # Output:
  #     thetas.sampled: p by p by nmc arrays of saved posterior samples of
  #     precision matrix
  #     lambdas.sampled: p*(p-1)/2 by nmc vector of saved samples of lambda
  #     squared (local tuning parameter)
  #     taus.sampled: 1 by nmc vector of saved samples of tau squared (global
  #     tuning parameter)
  p = nrow(S)
  omega_save = array(0,c(p,p,nmc))
  lambda_sq_save = matrix(0,p*(p-1)/2,nmc)
  tau = rep(0,nmc)
  tau_sq_save = rep(0,nmc)
  
  ind_all = matrix(0,p-1,p)
  for (i in 1:p){
    if (i==1) ind = matrix(2:p,ncol=1)
    else if (i==p) ind = matrix(1:(p-1),ncol=1)
    else ind = matrix(c(1:(i-1),(i+1):p),ncol=1)
    
    ind_all[,i] = ind
  }
  
  # set initial values
  Omega = diag(1,p) 
  Sigma = diag(1,p)
  Nu = matrix(1,p,p)
  Lambda_sq = matrix(1,p,p)
  Nu[1:p,1:p] = 1
  tau_sq = 1
  xi = 1
  
  for (iter in 1:(burnin+nmc)){  
    
    ### sample Sigma and Omega=inv(Sigma)
    for (i in 1:p){
      ind = ind_all[,i]     
      Sigma_11 = Sigma[ind,ind]
      sigma_12 = Sigma[ind,i]
      sigma_22 = Sigma[i,i]
      s_21 = S[ind,i]
      s_22 = S[i,i]
      lambda_sq_12 = Lambda_sq[ind,i]
      nu_12 = Nu[ind,i]
      ## sample gamma and beta
      gamma = rgamma(1,shape=(n/2+1),scale=2/s_22)    # random gamma with shape=n/2+1, rate=s_22/2. Sample 1. 
      #gamma= 0.06757892
      inv_Omega_11 = Sigma_11 - sigma_12%*%t(sigma_12)/sigma_22
      inv_C = s_22*inv_Omega_11+diag(1/(lambda_sq_12*tau_sq))
      #if(i!=1) return(inv_C)
      inv_C_chol = chol(inv_C) 
      #cat(iter,' ', i, '\n')
      mu_i= tryCatch({
        -solve(inv_C,s_21)
      }, error = function(e) {
        return(list(inv_C,s_21))
      })
      #mu_i = -solve(inv_C,s_21)
      beta = mu_i+ solve(inv_C_chol,rnorm(p-1))
      #kkk=c(0.1156077, 1.7215140, 0.10147592)
      #beta = mu_i+ solve(inv_C_chol,kkk)
      omega_12 = beta
      omega_22 = gamma + t(beta)%*%inv_Omega_11%*%beta;
      ## sample lambda_sq and nu
      rate = omega_12^2/(2*tau_sq)+1/nu_12
      lambda_sq_12 = 1/rgamma(length(omega_12), shape=1, scale = 1/rate)   # random inv gamma with shape=1, rate=rate. 
      nu_12 = 1/rgamma(length(lambda_sq_12), shape=1, scale = 1/(1+1/lambda_sq_12))   # random inv gamma with shape=1, rate=1+1/lambda_sq_12
      #lambda_sq_12 = 1/c(4.059683, 1.923486, 1.999748)
      #nu_12 = 1/c(1.5817724, 0.4459362, 0.9780721)
      ## update Omega, Sigma, Lambda_sq, Nu
      Omega[i,ind] = omega_12
      Omega[ind,i] = omega_12
      Omega[i,i] = omega_22
      temp = inv_Omega_11%*%beta
      Sigma_11 = inv_Omega_11 + temp%*%t(temp)/gamma;
      sigma_12 = -temp/gamma
      sigma_22 = 1/gamma
      Sigma[ind,ind] = Sigma_11
      Sigma[i,i] = sigma_22
      Sigma[i,ind] = sigma_12
      Sigma[ind,i] = sigma_12
      Lambda_sq[i,ind] = lambda_sq_12
      Lambda_sq[ind,i] = lambda_sq_12
      Nu[i,ind] = nu_12
      Nu[ind,i] = nu_12
    }
    # sample tau_sq and xi
    omega_vector =  Omega[lower.tri(Omega)] # Only elements in Omega below the diagonal
    lambda_sq_vector = Lambda_sq[lower.tri(Lambda_sq)]
    rate = 1/xi + sum(omega_vector^2/(2*lambda_sq_vector))
    tau_sq = 1/rgamma(1, shape=(p*(p-1)/2+1)/2, scale = 1/rate) # inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate. Only sample one. 
    xi = 1/rgamma(1, shape=1, scale = 1/(1+1/tau_sq))   # inv gamma w/ shape=1, rate=1+1/tau_sq
    
    # save Omega, lambda_sq, tau_sq
    if (iter > burnin){         
      omega_save[,,(iter-burnin)] = Omega
      lambda_sq_save[,(iter-burnin)] = lambda_sq_vector
      tau_sq_save[(iter-burnin)] = tau_sq
    }
  }
  return(list(thetas.sampled = omega_save,lambdas.sampled = lambda_sq_save,
              taus.samples = tau_sq_save))
}