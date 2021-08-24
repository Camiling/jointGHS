devtools::install_github("richardli/EMGS", subdir = "EMGS")
library(EMGS)
library(fastGHS)
library(tailoredGlasso) # utils, download from https://github.com/Camiling/tailoredGlasso
library(huge)
library(network)
library(GGally)

# Compare ECM GHS to EMGS (ECM for SSL)

# Not able to download package

set.seed(2020)
g <- huge.generator(n=50,d=100,graph = 'scale-free')
X <- scale(g$data)

# Use ECMGHS
res <- fastGHS(X, epsilon = 1e-7, savepath = T)
# select the 2% edges with the largest partial correlation to get the correct sparsity
theta.est <- cov2cor(res$theta)
theta.est.off.diag <- theta.est
diag(theta.est.off.diag) <- NA
theta.est[which(abs(theta.est) < quantile(abs(theta.est.off.diag), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag), 0.98,na.rm = T)

# Use EMGS
v0s <- seq(0.001, 0.5, len = 40)
v1 <- 100
lambda <- 1
a <- 2
b <- nrow(X)
res.emgs = EMGS(X, v0s, v1, lambda, a, b, epsilon = 1e-5, a_tau = 2, b_tau = 1)


