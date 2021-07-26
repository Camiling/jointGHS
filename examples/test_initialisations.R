
# Test different initialisations for fastGHS

# Not the most recent, as we do not ECMGHS select the sparsity itself here. See fixing_tau

# Initialise with the graphical lasso ---------------------

library(fastGHS)
library(tailoredGlasso) # utils, download from https://github.com/Camiling/tailoredGlasso
library(huge)

set.seed(2020)
g <- huge.generator(n=50,d=100,graph = 'scale-free')
X <- scale(g$data)
g$sparsity
# Use the graphical lasso
huge.g = huge(X,method='glasso',cov.output = T, lambda = 0.35737)
huge.g$sparsity
# Perform ECM for GHS with glasso graph as prior
res <- fastGHS(X, theta=huge.g$icov[[1]],epsilon = 1e-3, save_Q = T)
theta.est <- cov2cor(res$theta)
theta.est.off.diag <- theta.est
diag(theta.est.off.diag) <- NA
theta.est[which(abs(theta.est) < quantile(abs(theta.est.off.diag), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est!=0)
# 0.02

# ECHGHS
res.old <- fastGHS(X, epsilon = 1e-7)
theta.est.old <- cov2cor(res.old$theta)
theta.est.off.diag.old  <- theta.est.old 
diag(theta.est.off.diag.old ) <- NA
theta.est.old [which(abs(theta.est.old ) < quantile(abs(theta.est.off.diag.old ), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.old !=0)
#0.02
# precision of estimate
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est!=0)
# 0.3333333

# Mean difference between estimates
mean(abs(res$theta-res.old$theta))
# 6.838905e-14

# Precision wrt old estimate
tailoredGlasso::precision(theta.est!=0,theta.est.old!=0)
# 1

# Very similar!

# Initialise with glasso cov and icov ---------------------

# Perform ECM for GHS with glasso graph as prior
res2 <- fastGHS(X, theta=huge.g$icov[[1]], sigma=huge.g$cov[[1]],epsilon = 1e-7, save_Q = T)
# select the 2% edges with the largest partial correlation to get the correct sparsity
theta.est2 <- cov2cor(res2$theta)
theta.est.off.diag2 <- theta.est2
diag(theta.est.off.diag2) <- NA
theta.est2[which(abs(theta.est2) < quantile(abs(theta.est.off.diag2), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est2!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est2!=0)
# 0.3333333

mean(abs(res$theta-res2$theta))
# 6.526208e-19
tailoredGlasso::precision(theta.est!=0,theta.est2!=0)
# 1

# Similar, but not identical

# Initialise with glasso cov and icov, and different tau and lambda ---------------------

# Perform ECM for GHS with glasso graph as prior
lambda_sq_prior = matrix(rep(0.001,ncol(X)^2), ncol=ncol(X))
lambda_sq_prior[which(huge.g$icov[[1]]!=0, arr.ind=T)] = 5
res3 <- fastGHS(X, theta=huge.g$icov[[1]], sigma=huge.g$cov[[1]],tau_sq = 3,Lambda_sq=lambda_sq_prior ,epsilon = 1e-7)
theta.est3 <- cov2cor(res3$theta)
theta.est.off.diag3 <- theta.est3
diag(theta.est.off.diag3) <- NA
theta.est3[which(abs(theta.est3) < quantile(abs(theta.est.off.diag3), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est3!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est3!=0)
# 0.3333333
# same precision

mean(abs(res$theta-res3$theta))
# 9.855142e-14
tailoredGlasso::precision(theta.est!=0,theta.est3!=0)
# 0.969697

# Initialise with glasso cov and icov, and different tau but not lambda ---------------------

# Perform ECM for GHS with glasso graph as prior
res4 <- fastGHS(X, theta=huge.g$icov[[1]], sigma=huge.g$cov[[1]],tau_sq = 5,epsilon = 1e-7)
theta.est4 <- cov2cor(res4$theta)
theta.est.off.diag4 <- theta.est4
diag(theta.est.off.diag4) <- NA
theta.est4[which(abs(theta.est4) < quantile(abs(theta.est.off.diag4), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est4!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est4!=0)
# 0.3333333
# same precision

mean(abs(res$theta-res4$theta))
# 1.152475e-13
tailoredGlasso::precision(theta.est!=0,theta.est4!=0)
# 1


# Initialise with glasso icov, inverted glasso icov for sigma and different tau at start ---------------------

# Perform ECM for GHS with glasso graph as prior
res5 <- fastGHS(X, theta=huge.g$icov[[1]], sigma=solve(huge.g$icov[[1]]),tau_sq = 5,epsilon = 1e-7)
theta.est5 <- cov2cor(res5$theta)
theta.est.off.diag5 <- theta.est5
diag(theta.est.off.diag5) <- NA
theta.est5[which(abs(theta.est5) < quantile(abs(theta.est.off.diag5), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est5!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est5!=0)
# 0.3333333
# same precision

mean(abs(res$theta-res5$theta))
# 1.152475e-13
tailoredGlasso::precision(theta.est!=0,theta.est5!=0)
# 1

# Initialise with completely different values ---------------------

theta_init = matrix(runif(ncol(X)^2,0,1), ncol=ncol(X))
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res6 <- fastGHS(X, theta=theta_init, sigma=solve(theta_init),tau_sq = 1,epsilon = 1e-7, save_Q = T)
theta.est6 <- cov2cor(res6$theta)
theta.est.off.diag6 <- theta.est6
diag(theta.est.off.diag6) <- NA
theta.est6[which(abs(theta.est6) < quantile(abs(theta.est.off.diag6), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est6!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est6!=0)
# 0.2424242
# same precision

mean(abs(res$theta-res6$theta))
# 9.666188e-14
tailoredGlasso::precision(theta.est!=0,theta.est6!=0)
# 0.7070707



# Test with fixed Tau_sq --------------------------------

# Done by commenting out update for Tau
res.t <- fastGHS(X,tau_sq = 0.1 ,epsilon = 1e-10, save_Q = T)
theta.est.t <- cov2cor(res.t$theta)
theta.est.off.diag.t <- theta.est.t
diag(theta.est.off.diag.t) <- NA
theta.est.t[which(abs(theta.est.t) < quantile(abs(theta.est.off.diag.t), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.t), 0.98,na.rm = T) # the quantile
# 0.179965
tailoredGlasso::sparsity(theta.est.t!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est.t!=0)
# 0.2424242
theta.est.t[1:5,1:5]
#     [,1]      [,2]      [,3] [,4] [,5]
#[1,]    1 0.0000000 0.0000000    0    0
#[2,]    0 1.0000000 0.3794336    0    0
#[3,]    0 0.3794336 1.0000000    0    0
#[4,]    0 0.0000000 0.0000000    1    0
#[5,]    0 0.0000000 0.0000000    0    1

# Now elements are not as small! Could indicate an issue with its selection. Now try the others

# Test with fixed Lambda_sq --------------------------------

# Done by commenting out update for Lambda_sq
lambda_sq_prior = matrix(rep(0.1,ncol(X)^2), ncol=ncol(X))
lambda_sq_prior[which(huge.g$icov[[1]]!=0, arr.ind=T)] = 2
res.l <- fastGHS(X,Lambda_sq = lambda_sq_prior, epsilon = 1e-10, save_Q = T)
theta.est.l <- cov2cor(res.l$theta)
theta.est.off.diag.l <- theta.est.l
diag(theta.est.off.diag.l) <- NA
theta.est.l[which(abs(theta.est.l) < quantile(abs(theta.est.off.diag.l), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.l), 0.98,na.rm = T) # The quantile
# 8.012175e-17 
tailoredGlasso::sparsity(theta.est.l!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est.l!=0)
# 0.3333333
theta.est.l[1:5,1:5]
#             [,1]         [,2]        [,3] [,4]         [,5]
#[1,] 1.000000e+00 1.210606e-15 0.00000e+00    0 0.000000e+00
#[2,] 1.210606e-15 1.000000e+00 1.92322e-15    0 1.232505e-15
#[3,] 0.000000e+00 1.923220e-15 1.00000e+00    0 0.000000e+00
#[4,] 0.000000e+00 0.000000e+00 0.00000e+00    1 0.000000e+00
#[5,] 0.000000e+00 1.232505e-15 0.00000e+00    0 1.000000e+00

# Elements shrink a lot anyways!

res6$tau_sq
# 8.354465e-13


# Try with a more clever initial value of tau_sq, theta and lambda ---------------------

# Larger and informative values of theta, to make initial value fo tau larger. And more informative lambda
theta_init = as.matrix(Matrix::nearPD(2*(huge.g$icov[[1]]!=0))$mat)
lambda_init = matrix(rep(0.05, ncol(X)^2), ncol=ncol(X))
lambda_init[which(huge.g$icov[[1]]!=0, arr.ind=T)] = 1
res7 <- fastGHS(X, theta=theta_init, sigma=solve(theta_init),tau_sq = 1,Lambda_sq = lambda_init, epsilon = 1e-7, save_Q = T)
theta.est7 <- cov2cor(res7$theta)
theta.est.off.diag7 <- theta.est7
diag(theta.est.off.diag7) <- NA
theta.est7[which(abs(theta.est7) < quantile(abs(theta.est.off.diag7), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag7), 0.98,na.rm = T) # The quantile
tailoredGlasso::sparsity(theta.est7!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est7!=0)
# 0.3333333
# same precision


mean(abs(res$theta-res7$theta))
# 9.903212e-14
tailoredGlasso::precision(theta.est!=0,theta.est7!=0)
# 0.969697

# The first update for tau (if lambda was not updated first): 
#(2*(sum(theta_init^2/lambda_init)-sum(diag(theta_init^2/lambda_init)))/2+4/2)/(100*99+6)


# Inform initial values by truth ---------------------

# Larger and informative values of theta, to make initial value fo tau larger. And more informative lambda
theta_init = as.matrix(Matrix::nearPD(as.matrix(abs(g$omega)>1e-7))$mat)
lambda_init = matrix(rep(0.05, ncol(X)^2), ncol=ncol(X))
lambda_init[which(as.matrix(abs(g$omega)>1e-7), arr.ind=T)] = 1
res8 <- fastGHS(X, theta=theta_init, sigma=solve(theta_init),tau_sq = 1,Lambda_sq = lambda_init, epsilon = 1e-7, save_Q = T)
theta.est8 <- cov2cor(res8$theta)
theta.est.off.diag8 <- theta.est8
diag(theta.est.off.diag8) <- NA
theta.est8[which(abs(theta.est8) < quantile(abs(theta.est.off.diag8), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag8), 0.98, na.rm = T) # The quantile
# 1.647206e-15 
tailoredGlasso::sparsity(theta.est8!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est8!=0)
# 0.959596
# same precision


mean(abs(res$theta-res8$theta))
# 1.00294e-13
tailoredGlasso::precision(theta.est!=0,theta.est8!=0)
# 0.3737374

theta.est8[1:5,1:5]
#[,1]         [,2]         [,3]         [,4]         [,5]
#[1,] 1.000000e+00 2.614752e-15 0.000000e+00 0.000000e+00 0.000000e+00
#[2,] 2.614752e-15 1.000000e+00 2.145038e-14 3.622705e-16 2.441853e-15
#[3,] 0.000000e+00 2.145038e-14 1.000000e+00 0.000000e+00 0.000000e+00
#[4,] 0.000000e+00 3.622705e-16 0.000000e+00 1.000000e+00 0.000000e+00
#[5,] 0.000000e+00 2.441853e-15 0.000000e+00 0.000000e+00 1.000000e+00

# Prec mat elements are still very small, but the precision is very high!!



# Compare to standard GHS --------------------------------------


source('~/Documents/Cambridge/PhD/Graph reconstruction/GHS/GHS.R')

# GENERATE GRAPH 4: n=100, p=50, larger partial correlations (0.2)
n.4=100
p.4=50
set.seed(12345)
data.sf.4 = huge.generator(n=n.4, d=p.4,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.4 = data.sf.4$theta # True adjacency matrix
theta.true.4 = data.sf.4$omega # The precision matrix
theta.true.4[which(theta.true.4<10e-5,arr.ind=T)]=0  
g.sf.4=graph.adjacency(data.sf.4$theta,mode="undirected",diag=F) # true igraph object
x.sf.4 = data.sf.4$data # Observed attributes. nxp matrix.
x.sf.scaled.4= scale(x.sf.4) # Scale columns/variables.
s.sf.scaled.4 = cov(x.sf.scaled.4) # Empirical covariance matrix
data.sf.4$sparsity # True sparsity: 0.04

# Ordinary GHS
set.seed(123)
ghs.res.4 = GHS(t(x.sf.4)%*%x.sf.4*n.4,n.4,burnin=100,nmc=1000)
hist(ghs.res.4$taus.samples) # Small tau, around size 1e-8
theta.est.ghs = cov2cor(apply(ghs.res.4$thetas.sampled, c(1,2), mean))
theta.est.off.diag.ghs <- theta.est.ghs
diag(theta.est.off.diag.ghs) <- NA
theta.est.ghs[which(abs(theta.est.ghs) < quantile(abs(theta.est.off.diag.ghs), 0.96,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.ghs), 0.96,na.rm = T)
# 0.05488769 
tailoredGlasso::sparsity(theta.est.ghs!=0)
# 0.04
tailoredGlasso::precision(theta.true.4!=0,theta.est.ghs!=0)
# 0.7142857

theta.est.ghs[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1398641 0.0000000 0.0000000 0.0000000
#[2,] 0.1398641 1.0000000 0.2282595 0.0000000 0.0000000
#[3,] 0.0000000 0.2282595 1.0000000 0.1833118 0.2420606
#[4,] 0.0000000 0.0000000 0.1833118 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2420606 0.0000000 1.0000000

# Closer look at distribution of taus
hist(ghs.res.4$taus.samples,breaks=100)

# Closer look at distribution of lambdas. Diagonal elements are not included here
hist(ghs.res.4$lambdas.sampled[which(ghs.res.4$lambdas.sampled[,1000]<1),1000],breaks=100)
# look at the smaller ones
hist(ghs.res.4$lambdas.sampled[which(ghs.res.4$lambdas.sampled[,1000]<1e-3),1000],breaks=50)
summary(ghs.res.4$lambdas.sampled[,1000])
sum(ghs.res.4$lambdas.sampled>10)
# 238441

# In GHS, many lambdas are very large!

# ECM GHS

res9 <- fastGHS(x.sf.4, tau_sq = 1, epsilon = 1e-7, save_Q = T)
theta.est9 <- cov2cor(res9$theta)
theta.est.off.diag9 <- theta.est9
diag(theta.est.off.diag9) <- NA
theta.est9[which(abs(theta.est9) < quantile(abs(theta.est.off.diag9), 0.96,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag9), 0.96, na.rm = T) # The quantile
# 1.092516e-22 
tailoredGlasso::sparsity(theta.est9!=0)
# 0.04
tailoredGlasso::precision(theta.true.4!=0, theta.est9!=0)
# 0.6530612

theta.est9[1:5,1:5]
#[,1]         [,2]         [,3]         [,4]         [,5]
#[1,] 1.000000e+00 1.298145e-19 0.000000e+00 0.000000e+00 0.000000e+00
#[2,] 1.298145e-19 1.000000e+00 3.035156e-13 0.000000e+00 0.000000e+00
#[3,] 0.000000e+00 3.035156e-13 1.000000e+00 2.214765e-15 3.903557e-14
#[4,] 0.000000e+00 0.000000e+00 2.214765e-15 1.000000e+00 0.000000e+00
#[5,] 0.000000e+00 0.000000e+00 3.903557e-14 0.000000e+00 1.000000e+00

# Worse precision for ECM than GHS.But GHS rarely applicable for even medium networks. 

res9$tau_sq

diag(res9$Lambda_sq)=0
hist(res9$Lambda_sq[which(res9$Lambda_sq<0.02 & res9$Lambda_sq!=0)],breaks=50)
max(res9$Lambda_sq)
min(res9$Lambda_sq[which(res9$Lambda_sq!=0)])


# Test with fixed Tau_sq on this same data --------------------------------

res.10 <- fastGHS(x.sf.4,tau_sq = 0.1 ,epsilon = 1e-5, save_Q = T, fix_tau=TRUE)
theta.est.10 <- cov2cor(res.10$theta)
theta.est.off.diag.10 <- theta.est.10
diag(theta.est.off.diag.10) <- NA
theta.est.10[which(abs(theta.est.10) < quantile(abs(theta.est.off.diag.10), 0.96,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.10), 0.96,na.rm = T)
# 98% 
# 8.5503e-24
tailoredGlasso::sparsity(theta.est.10!=0)
# 0.04
tailoredGlasso::precision(as.matrix(theta.true.4!=0), theta.est.10!=0)
# 0.6530612
theta.est.10[1:5,1:5]
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2026686 0.0000000 0.0000000 0.0000000
#[2,] 0.2026686 1.0000000 0.2464950 0.0000000 0.0000000
#[3,] 0.0000000 0.2464950 1.0000000 0.2176918 0.2528399
#[4,] 0.0000000 0.0000000 0.2176918 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2528399 0.0000000 1.0000000

# Same precision

# A bit better as we do not get extreme underflow, at same precision. 

# Test with fixed Tau_sq again --------------------------------

res.t2 <- fastGHS(X,tau_sq = 0.001 ,epsilon = 1e-7, save_Q = T, fix_tau=TRUE)
theta.est.t2 <- cov2cor(res.t2$theta)
theta.est.off.diag.t2 <- theta.est.t2
diag(theta.est.off.diag.t2) <- NA
theta.est.t2[which(abs(theta.est.t2) < quantile(abs(theta.est.off.diag.t2), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.t2), 0.98,na.rm = T)
# 98% 
# 2.344093e-08 
tailoredGlasso::sparsity(theta.est.t2!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est.t2!=0)
# 0.3333333
theta.est.t2[1:5,1:5]
#[,1]         [,2]       [,3] [,4]         [,5]
#[1,] 1.000000e+00 2.561135e-08 0.0000e+00    0 0.000000e+00
#[2,] 2.561135e-08 1.000000e+00 4.6127e-08    0 2.615184e-08
#[3,] 0.000000e+00 4.612700e-08 1.0000e+00    0 0.000000e+00
#[4,] 0.000000e+00 0.000000e+00 0.0000e+00    1 0.000000e+00
#[5,] 0.000000e+00 2.615184e-08 0.0000e+00    0 1.000000e+00


# A bit better as we do not get extreme underflow, as same precision. 

# Note that in the above examples, we force too many edges to be included: the method decides the sparsity itself, and so a simple thresholding rule is sufficient
