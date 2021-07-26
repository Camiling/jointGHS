library(fastGHS)
library(huge)
library(glasso)
library(igraph)
source('GHS.R')

# Compare results to those from the GHS-like prior (Sagar et al., 2021) (preprint). 

# EXAMPLE 1 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=50, larger partial correlations (0.229)

n.test=100
p.test=50
set.seed(12345)
data.sf.test = huge::huge.generator(n=n.test, d=p.test,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.test = data.sf.test$theta # True adjacency matrix
theta.true.test = data.sf.test$omega # The precision matrix
theta.true.test[which(theta.true.test<10e-5,arr.ind=T)]=0  
g.sf.test=graph.adjacency(data.sf.test$theta,mode="undirected",diag=F) # true igraph object
x.sf.test = data.sf.test$data # Observed attributes. nxp matrix.
data.sf.test$sparsity # True sparsity: 0.04
# Look at precision matrix (partial correlations)
cov2cor(theta.true.test[1:5,1:5])

# ECMGHS with fixed tau

res.test <- fastGHS(x.sf.test,tau_sq = 0.1,epsilon = 1e-3, fix_tau=TRUE)
theta.est.test <- cov2cor(res.test$theta)
theta.est.test[which(abs(theta.est.test) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.test!=0)
# 0.01632653
tailoredGlasso::precision(as.matrix(theta.true.test!=0), theta.est.test!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.test!=0), theta.est.test!=0)
# 0.4081633
theta.est.test[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1965087 0.0000000 0.0000000 0.0000000
#[2,] 0.1965087 1.0000000 0.2468123 0.0000000 0.0000000
#[3,] 0.0000000 0.2468123 1.0000000 0.2176918 0.2528399
#[4,] 0.0000000 0.0000000 0.2176918 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2528399 0.0000000 1.0000000

# Test ECMGHS for several tau values

taus = seq(1e-6,0.5,length.out = 300)
res.test.ghs.list <- lapply(taus,FUN = function(l) fastGHS(x.sf.test,tau_sq = l,epsilon = 1e-3, fix_tau = T))
plot(taus, unlist(lapply(res.test.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))))
# As we see, the sparsity plateaus
plot(taus, unlist(lapply(res.test.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.test!=0),abs(s$theta)>1e-5))))


# GHS-like prior

set.seed(123)
res.test.ghslike <- fastGHS(x.sf.test,tau_sq = 0.001,epsilon = 1e-5, GHS_like=TRUE) # Now a is tau_sq
theta.est.test.ghslike <- cov2cor(res.test.ghslike$theta)
theta.est.test.ghslike[which(abs(theta.est.test.ghslike) < quantile(abs(theta.est.test.ghslike)[upper.tri(theta.est.test.ghslike)], 0.96), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.test.ghslike!=0)
# 0.04
tailoredGlasso::precision(as.matrix(theta.true.test!=0), theta.est.test.ghslike!=0)
# 0.3061224
tailoredGlasso::recall(as.matrix(theta.true.test!=0), theta.est.test.ghslike!=0)
# 0.3061224
theta.est.test.ghslike[1:5,1:5]
#[,1]          [,2]         [,3]          [,4]          [,5]
#[1,]    1  0.000000e+00 0.000000e+00 0.000000e+00  0.000000e+00
#[2,]    0  1.000000e+00 6.710463e-15 0.000000e+00 -4.510317e-16
#[3,]    0  6.710463e-15 1.000000e+00 1.808246e-15  8.812638e-16
#[4,]    0  0.000000e+00 1.808246e-15 1.000000e+00  0.000000e+00
#[5,]    0 -4.510317e-16 8.812638e-16 0.000000e+00  1.000000e+00

# Now all elements shrink too much.... As in our original problem.

# Also, there is no clear distinguishment between almost-zero and non-zero elements. 

# Test for several tau values
taus = seq(1e-5,0.2,length.out = 200)
res.test.ghslike.list <- lapply(taus,FUN = function(l) fastGHS(x.sf.test,tau_sq = l,epsilon = 1e-5, GHS_like=TRUE))
plot(taus, unlist(lapply(res.test.ghslike.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))))

# As we see, we get no elements larger than 1e-5 for any choice of tau.


