library(jointGHS)
library(huge)
library(glasso)
library(igraph)
library(fastGHS)
source('examples/help_functions.R')

## Simple test

# EXAMPLE 1: two datasets from the same distribution ------------------------------------------------------------------

# n=100, p=50, larger partial correlations (0.229)
n.1=100
p.1=50
set.seed(12345)
data.sf.1= huge::huge.generator(n=n.1, d=p.1,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.1 = data.sf.1$theta # True adjacency matrix
theta.true.1 = data.sf.1$omega # The precision matrix
theta.true.1[which(theta.true.1<10e-5,arr.ind=T)]=0  
g.sf.1=graph.adjacency(data.sf.1$theta,mode="undirected",diag=F) # true igraph object
x.sf.1 = data.sf.1$data # Observed attributes. nxp matrix.
x.sf.scaled.1= scale(x.sf.1) # Scale columns/variables.
s.sf.scaled.1 = cov(x.sf.scaled.1) # Empirical covariance matrix
data.sf.1$sparsity # True sparsity: 0.04
cov2cor(theta.true.1)[1:5,1:5]          
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2290578 0.0000000 0.0000000 0.0000000
#[2,] 0.2290578 1.0000000 0.2290578 0.0000000 0.0000000
#[3,] 0.0000000 0.2290578 1.0000000 0.2290578 0.2290578
#[4,] 0.0000000 0.0000000 0.2290578 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2290578 0.0000000 1.0000000

# Generate second data set from the same distribution
n.1.2 =  100
set.seed(123)
x2.sf.1 = mvtnorm::rmvnorm(n.1.2, sigma = data.sf.1$sigma)
x2.sf.1.scaled = scale(x2.sf.1)

# Use jointGHS on two identical data sets
res.joint.1 = jointGHS::jointGHS(list(x.sf.scaled.1, x2.sf.1.scaled), tau_sq=c(10, 10),epsilon = 1e-3, fix_tau=TRUE)

theta1.est.1 <- cov2cor(res.joint.1$theta[[1]])
theta2.est.1 <- cov2cor(res.joint.1$theta[[2]])
theta1.est.1[which(abs(theta1.est.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta1.est.1!=0)
# 0.0122449
theta2.est.1[which(abs(theta2.est.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta2.est.1!=0)
# 0.0122449
tailoredGlasso::precision(as.matrix(theta.true.1!=0), theta1.est.1!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.1!=0), theta2.est.1!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.1!=0), theta1.est.1!=0)
# 0.3061224
tailoredGlasso::recall(as.matrix(theta.true.1!=0), theta1.est.1!=0)
# 0.3061224
theta1.est.1[1:5,1:5]
#[,1]      [,2]      [,3] [,4]      [,5]
#[1,]    1 0.0000000 0.0000000    0 0.0000000
#[2,]    0 1.0000000 0.2983643    0 0.0000000
#[3,]    0 0.2983643 1.0000000    0 0.2903136
#[4,]    0 0.0000000 0.0000000    1 0.0000000
#[5,]    0 0.0000000 0.2903136    0 1.0000000


# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.1 <- fastGHS(x.sf.scaled.1,tau_sq = 0.024,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.1 <- cov2cor(res.ecmghs.1$theta)
theta.est.ecmghs.1[which(abs(theta.est.ecmghs.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.1!=0)
# 0.0122449
tailoredGlasso::precision(as.matrix(theta.true.1!=0), theta.est.ecmghs.1!=0)
# 0.9333333
tailoredGlasso::recall(as.matrix(theta.true.1!=0), theta.est.ecmghs.1!=0)
# 0.2857143

set.seed(22)
res.ecmghs.1.2 <- fastGHS(x2.sf.1.scaled,tau_sq = 0.022,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.1.2 <- cov2cor(res.ecmghs.1.2$theta)
theta.est.ecmghs.1.2[which(abs(theta.est.ecmghs.1.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.1.2!=0)
# 0.0122449
tailoredGlasso::precision(as.matrix(theta.true.1!=0), theta.est.ecmghs.1.2!=0)
# 0.9333333
tailoredGlasso::recall(as.matrix(theta.true.1!=0), theta.est.ecmghs.1.2!=0)
# 0.2857143

# JointGHS outperforms the single-network version in terms of precision and recall. 

# EXAMPLE 2: two datasets from two completely unrelated distributions ------------------------------------------------------------------

# n=100, p=50, larger partial correlations (0.229)
n.2.1=100
p.2.1=50
set.seed(12345)
data.sf.2.1= huge::huge.generator(n=n.2.1, d=p.2.1,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.2.1 = data.sf.2.1$theta # True adjacency matrix
theta.true.2.1 = data.sf.2.1$omega # The precision matrix
theta.true.2.1[which(theta.true.2.1<10e-5,arr.ind=T)]=0  
g.sf.2.1=graph.adjacency(data.sf.2.1$theta,mode="undirected",diag=F) # true igraph object
x.sf.2.1 = data.sf.2.1$data # Observed attributes. nxp matrix.
x.sf.scaled.2.1= scale(x.sf.2.1) # Scale columns/variables.
s.sf.scaled.2.1 = cov(x.sf.scaled.2.1) # Empirical covariance matrix
data.sf.2.1$sparsity # True sparsity: 0.04

# Generate second data set 
# n=200, p=50
n.2.2=200
p.2.2=50
set.seed(123456)
data.sf.2.2= huge::huge.generator(n=n.2.2, d=p.2.2,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.2.2 = data.sf.2.2$theta # True adjacency matrix
theta.true.2.2 = data.sf.2.2$omega # The precision matrix
theta.true.2.2[which(theta.true.2.2<10e-5,arr.ind=T)]=0  
g.sf.2.2=graph.adjacency(data.sf.2.2$theta,mode="undirected",diag=F) # true igraph object
x.sf.2.2 = data.sf.2.2$data # Observed attributes. nxp matrix.
x.sf.scaled.2.2= scale(x.sf.2.2) # Scale columns/variables.
s.sf.scaled.2.2 = cov(x.sf.scaled.2.2) # Empirical covariance matrix
data.sf.2.2$sparsity # True sparsity: 0.04

# Use jointGHS on two unrelated data sets
res.joint.2 = jointGHS::jointGHS(list(x.sf.scaled.2.1, x.sf.scaled.2.2), tau_sq=c(10, 10),epsilon = 1e-3, fix_tau=TRUE)

theta1.est.2 <- cov2cor(res.joint.2$theta[[1]])
theta2.est.2 <- cov2cor(res.joint.2$theta[[2]])
theta1.est.2[which(abs(theta1.est.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta1.est.2!=0)
# 0.006530612
theta2.est.2[which(abs(theta2.est.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta2.est.2!=0)
# 0.005714286
tailoredGlasso::precision(as.matrix(theta.true.2.1!=0), theta1.est.2!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.2.2!=0), theta2.est.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.2.1!=0), theta1.est.2!=0)
# 0.1632653
tailoredGlasso::recall(as.matrix(theta.true.2.2!=0), theta2.est.2!=0)
# 0.1428571

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.2.1 <- fastGHS(x.sf.scaled.2.1,tau_sq = 0.008,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.2.1 <- cov2cor(res.ecmghs.2.1$theta)
theta.est.ecmghs.2.1[which(abs(theta.est.ecmghs.2.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.2.1!=0)
# 0.006530612
tailoredGlasso::precision(as.matrix(theta.true.2.1!=0), theta.est.ecmghs.2.1!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.2.1!=0), theta.est.ecmghs.2.1!=0)
# 0.1632653

set.seed(22)
res.ecmghs.2.2 <- fastGHS(x.sf.scaled.2.2,tau_sq = 0.003,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.2.2 <- cov2cor(res.ecmghs.2.2$theta)
theta.est.ecmghs.2.2[which(abs(theta.est.ecmghs.2.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.2.2!=0)
# 0.005714286
tailoredGlasso::precision(as.matrix(theta.true.2.2!=0), theta.est.ecmghs.2.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.2.2!=0), theta.est.ecmghs.2.2!=0)
# 0.1428571

# For competely unrelated networks, jointGHS gives the same results as the ordinary ECM GHS.

# EXAMPLE 3: two datasets from related distributions ------------------------------------------------------------------

# 10% edge disagreement in the underlying edges

# First network: n=100, p=50, larger partial correlations (0.229)
n.3.1=100
p.3.1=50
set.seed(12345)
data.sf.3.1= huge::huge.generator(n=n.3.1, d=p.3.1,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.3.1 = data.sf.3.1$theta # True adjacency matrix
theta.true.3.1 = data.sf.3.1$omega # The precision matrix
theta.true.3.1[which(theta.true.3.1<10e-5,arr.ind=T)]=0  
g.sf.3.1=graph.adjacency(data.sf.3.1$theta,mode="undirected",diag=F) # true igraph object
x.sf.3.1 = data.sf.3.1$data # Observed attributes. nxp matrix.
x.sf.scaled.3.1= scale(x.sf.3.1) # Scale columns/variables.
s.sf.scaled.3.1 = cov(x.sf.scaled.3.1) # Empirical covariance matrix
data.sf.3.1$sparsity # True sparsity: 0.04

# Generate second data set with same sparsity
# n=200, p=50, 10% edge disagreement
n.3.2=200
p.3.2=50
set.seed(123456)
graph.3.2 = mutate.graph(data.sf.3.1, 0.1)
theta.true.3.2 = graph.3.2$prec.mat
tailoredGlasso::sparsity(theta.true.3.2!=0)
# 0.04
x.sf.scaled.3.2 = scale(mvtnorm::rmvnorm(n.3.2, sigma = solve(graph.3.2$prec.mat)))

# Use jointGHS on two unrelated data sets
res.joint.3 = jointGHS::jointGHS(list(x.sf.scaled.3.1, x.sf.scaled.3.2), tau_sq=c(10, 10),epsilon = 1e-3, fix_tau=TRUE)

theta1.est.3 <- cov2cor(res.joint.3$theta[[1]])
theta2.est.3 <- cov2cor(res.joint.3$theta[[2]])
theta1.est.3[which(abs(theta1.est.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta1.est.3!=0)
# 0.01061224
theta2.est.3[which(abs(theta2.est.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta2.est.3!=0)
# 0.01142857
tailoredGlasso::precision(as.matrix(theta.true.3.1!=0), theta1.est.3!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.3.2!=0), theta2.est.3!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.3.1!=0), theta1.est.3!=0)
# 0.2653061
tailoredGlasso::recall(as.matrix(theta.true.3.2!=0), theta2.est.3!=0)
# 0.2857143

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.3.1 <- fastGHS(x.sf.scaled.3.1,tau_sq = 0.017,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.3.1 <- cov2cor(res.ecmghs.3.1$theta)
theta.est.ecmghs.3.1[which(abs(theta.est.ecmghs.3.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.3.1!=0)
# 0.01061224
tailoredGlasso::precision(as.matrix(theta.true.3.1!=0), theta.est.ecmghs.3.1!=0)
# 0.9230769
tailoredGlasso::recall(as.matrix(theta.true.3.1!=0), theta.est.ecmghs.3.1!=0)
# 0.244898

set.seed(22)
res.ecmghs.3.2 <- fastGHS(x.sf.scaled.3.2,tau_sq = 0.01,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.3.2 <- cov2cor(res.ecmghs.3.2$theta)
theta.est.ecmghs.3.2[which(abs(theta.est.ecmghs.3.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.3.2!=0)
# 0.01142857
tailoredGlasso::precision(as.matrix(theta.true.3.2!=0), theta.est.ecmghs.3.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.3.2!=0), theta.est.ecmghs.3.2!=0)
# 0.2857143

# For related networks, jointGHS gives results similar to the ordinary ECM GHS, but slightly more accurate. 

# Slightly better precision and recall than the single-network version



res.joint.3.2 = jointGHS::jointGHS(list(x.sf.scaled.3.1, x.sf.scaled.3.2), tau_sq=c(1000, 1000),epsilon = 1e-3, fix_tau=TRUE)
res.joint.3.2$epsilon





