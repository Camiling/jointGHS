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

# ECMGHS

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


# GHS-like prior

res.test.ghslike <- fastGHS(x.sf.test,tau_sq = 0.001,epsilon = 1e-5, GHS_like=TRUE) # Now a is tau_sq
theta.est.test.ghslike <- cov2cor(res.test.ghslike$theta)
theta.est.test.ghslike[which(abs(theta.est.test.ghslike) < quantile(abs(theta.est.test.ghslike)[upper.tri(theta.est.test.ghslike)], 0.96), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.test.ghslike!=0)
# 0.04
tailoredGlasso::precision(as.matrix(theta.true.test!=0), theta.est.test.ghslike!=0)
# 0.3265306
tailoredGlasso::recall(as.matrix(theta.true.test!=0), theta.est.test.ghslike!=0)
# 0.3265306
theta.est.test.ghslike[1:5,1:5]
#[,1]          [,2]         [,3]          [,4]          [,5]
#[1,]    1  0.000000e+00 0.000000e+00  0.000000e+00  0.000000e+00
#[2,]    0  1.000000e+00 5.977752e-14 -1.846785e-18 -7.155238e-18
#[3,]    0  5.977752e-14 1.000000e+00  2.515150e-18  1.071527e-14
#[4,]    0 -1.846785e-18 2.515150e-18  1.000000e+00 -3.877681e-17
#[5,]    0 -7.155238e-18 1.071527e-14 -3.877681e-17  1.000000e+00

# Now all elements shrink too much.... As in our original problem.

# Also, there is no clear distinguishment between almost-zero and non-zero elements. 

# Test for several tau values
taus = seq(1e-5,0.2,length.out = 1000)
res.test.ghslike.list <- lapply(taus,FUN = function(l) fastGHS(x.sf.test,tau_sq = l,epsilon = 1e-5, GHS_like=TRUE))
plot(taus, unlist(lapply(res.test.ghslike.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))))

# As we see, we get no elements larger than 1e-5 for any choice of tau.


