
# Test different initialisations for fastGHS

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
res <- fastGHS(X, theta=huge.g$icov[[1]],epsilon = 1e-7, save_Q = T)
theta.est <- cov2cor(res$theta)
theta.est.off.diag <- theta.est
diag(theta.est.off.diag) <- NA
theta.est[which(abs(theta.est) < quantile(abs(theta.est.off.diag), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est!=0)
# 0.02

# Compare to identity matrix as prior
res.old <- fastGHS(X, epsilon = 1e-7)
theta.est.old <- cov2cor(res.old$theta)
theta.est.off.diag.old  <- theta.est.old 
diag(theta.est.off.diag.old ) <- NA
theta.est.old [which(abs(theta.est.old ) < quantile(abs(theta.est.off.diag.old ), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.old !=0)
0.02
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
# 1

# Test with fixed Tau_sq --------------------------------

res.t <- fastGHS(X,tau_sq = 0.1 ,epsilon = 1e-10, save_Q = T)
theta.est.t <- cov2cor(res.t$theta)
theta.est.off.diag.t <- theta.est.t
diag(theta.est.off.diag.t) <- NA
theta.est.t[which(abs(theta.est.t) < quantile(abs(theta.est.off.diag.t), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.t!=0)
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est.t!=0)
# 0.2424242
theta.est.t[1:5,1:5]
#[,1]      [,2]      [,3] [,4] [,5]
#[1,]    1 0.0000000 0.0000000    0    0
#[2,]    0 1.0000000 0.3794336    0    0
#[3,]    0 0.3794336 1.0000000    0    0
#[4,]    0 0.0000000 0.0000000    1    0
#[5,]    0 0.0000000 0.0000000    0    1

# Now elements are not as small! Could indicate an issue with its selection. Now try the others

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




