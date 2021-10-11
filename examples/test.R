library(jointGHS)
library(huge)
library(glasso)
library(igraph)
library(fastGHS)
library(ggplot2)
library(gridExtra)
source('examples/help_functions.R')

# Explore the performance of jointGHS compared to the single-network version in different settings

# Example 1: two datasets from the same distribution
# Example 2: two datasets from two completely unrelated distributions
# Example 3: two datasets from related distributions
# Example 4: four datasets from related distributions, more high-dimensional
# Example 5: ten datasets from related distributions, more high-dimensional

# EXAMPLE 1: two datasets from the same distribution ------------------------------------------------------------------

# n=100, p=50, larger partial correlations (0.229)
n.1=100
p.1=50
set.seed(12345)
data.sf.1 = huge::huge.generator(n=n.1, d=p.1,graph = 'scale-free',v=0.5,u=0.05) 
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

# 20% edge disagreement in the underlying edges

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
# n=200, p=50, 20% edge disagreement
n.3.2=200
p.3.2=50
set.seed(123456)
graph.3.2 = mutate.graph(data.sf.3.1, 0.2)
theta.true.3.2 = graph.3.2$prec.mat
tailoredGlasso::sparsity(theta.true.3.2!=0)
# 0.04
x.sf.scaled.3.2 = scale(mvtnorm::rmvnorm(n.3.2, sigma = solve(graph.3.2$prec.mat)))

# Use jointGHS on two slightly related data sets
res.joint.3 = jointGHS::jointGHS(list(x.sf.scaled.3.1, x.sf.scaled.3.2), tau_sq=c(1e-4, 1e-4),epsilon = 1e-3, fix_tau=TRUE)

theta1.est.3 <- cov2cor(res.joint.3$theta[[1]])
theta2.est.3 <- cov2cor(res.joint.3$theta[[2]])
theta1.est.3[which(abs(theta1.est.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta1.est.3!=0)
# 0.01142857
theta2.est.3[which(abs(theta2.est.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta2.est.3!=0)
# 0.01142857
tailoredGlasso::precision(as.matrix(theta.true.3.1!=0), theta1.est.3!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.3.2!=0), theta2.est.3!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.3.1!=0), theta1.est.3!=0)
# 0.2857143
tailoredGlasso::recall(as.matrix(theta.true.3.2!=0), theta2.est.3!=0)
# 0.2857143

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.3.1 <- fastGHS(x.sf.scaled.3.1,tau_sq = 0.019,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.3.1 <- cov2cor(res.ecmghs.3.1$theta)
theta.est.ecmghs.3.1[which(abs(theta.est.ecmghs.3.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.3.1!=0)
# 0.01142857
tailoredGlasso::precision(as.matrix(theta.true.3.1!=0), theta.est.ecmghs.3.1!=0)
# 0.9285714
tailoredGlasso::recall(as.matrix(theta.true.3.1!=0), theta.est.ecmghs.3.1!=0)
# 0.2653061

set.seed(22)
res.ecmghs.3.2 <- fastGHS(x.sf.scaled.3.2,tau_sq = 0.016,epsilon = 1e-3, fix_tau=TRUE)
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


# EXAMPLE 4: four datasets from related distributions, more high-dimensional ------------------------------------------------------------------

# 20% edge disagreement with the first network in the underlying edges
# Note: the 80% agreement is only with graph 1. I.e. networks 2, 3 and 4 are not as similar to each other. 

# First network: n=100, p=100
n.4.1=100
p.4.1=100
set.seed(12345)
data.sf.4.1= huge::huge.generator(n=n.4.1, d=p.4.1,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.4.1 = data.sf.4.1$theta # True adjacency matrix
theta.true.4.1 = data.sf.4.1$omega # The precision matrix
theta.true.4.1[which(theta.true.4.1<10e-5,arr.ind=T)]=0  
g.sf.4.1=graph.adjacency(data.sf.4.1$theta,mode="undirected",diag=F) # true igraph object
x.sf.4.1 = data.sf.4.1$data # Observed attributes. nxp matrix.
x.sf.scaled.4.1= scale(x.sf.4.1) # Scale columns/variables.
s.sf.scaled.4.1 = cov(x.sf.scaled.4.1) # Empirical covariance matrix
data.sf.4.1$sparsity # True sparsity: 0.02

# Generate second data set with same sparsity
# n=200, p=100, 20% edge disagreement
n.4.2=200
p.4.2=100
set.seed(123456)
graph.4.2 = mutate.graph(data.sf.4.1, 0.2)
theta.true.4.2 = graph.4.2$prec.mat
tailoredGlasso::sparsity(theta.true.4.2!=0)
# 0.02
x.sf.scaled.4.2 = scale(mvtnorm::rmvnorm(n.4.2, sigma = solve(graph.4.2$prec.mat)))

# Generate third data set with same sparsity
# n=150, p=100, 20% edge disagreement
n.4.3=150
p.4.3=100
set.seed(12345677)
graph.4.3 = mutate.graph(data.sf.4.1, 0.2)
theta.true.4.3 = graph.4.3$prec.mat
tailoredGlasso::sparsity(theta.true.4.3!=0)
# 0.02
x.sf.scaled.4.3 = scale(mvtnorm::rmvnorm(n.4.3, sigma = solve(graph.4.3$prec.mat)))

# Generate fourth data set with same sparsity
# n=250, p=100, 20% edge disagreement
n.4.4=250
p.4.4=100
set.seed(12345688)
graph.4.4 = mutate.graph(data.sf.4.1, 0.2)
theta.true.4.4 = graph.4.4$prec.mat
tailoredGlasso::sparsity(theta.true.4.4!=0)
# 0.02
x.sf.scaled.4.4 = scale(mvtnorm::rmvnorm(n.4.4, sigma = solve(graph.4.4$prec.mat)))

# How much do the networks agree priorly?
tailoredGlasso::precision(theta.true.4.1!=0,theta.true.4.2!=0)
# 0.8383838
tailoredGlasso::precision(theta.true.4.1!=0,theta.true.4.3!=0)
# 0.8181818
tailoredGlasso::precision(theta.true.4.1!=0,theta.true.4.4!=0)
# 0.8383838
tailoredGlasso::precision(theta.true.4.2!=0,theta.true.4.3!=0)
# 0.6868687
tailoredGlasso::precision(theta.true.4.2!=0,theta.true.4.4!=0)
# 0.6868687
tailoredGlasso::precision(theta.true.4.3!=0,theta.true.4.4!=0)
# 0.6767677

# So similarity is not symmetric.

# Use jointGHS on the four related data sets
set.seed(1234)
res.joint.4 = jointGHS::jointGHS(list(x.sf.scaled.4.1, x.sf.scaled.4.2, x.sf.scaled.4.3, x.sf.scaled.4.4), tau_sq=c(10, 10, 10, 10),epsilon = 1e-3, fix_tau=TRUE)

theta1.est.4 <- cov2cor(res.joint.4$theta[[1]])
theta2.est.4 <- cov2cor(res.joint.4$theta[[2]])
theta3.est.4 <- cov2cor(res.joint.4$theta[[3]])
theta4.est.4 <- cov2cor(res.joint.4$theta[[4]])
theta1.est.4[which(abs(theta1.est.4) < 1e-5, arr.ind = T)] = 0
theta2.est.4[which(abs(theta2.est.4) < 1e-5, arr.ind = T)] = 0
theta3.est.4[which(abs(theta3.est.4) < 1e-5, arr.ind = T)] = 0
theta4.est.4[which(abs(theta4.est.4) < 1e-5, arr.ind = T)] = 0

tailoredGlasso::sparsity(theta1.est.4!=0)
# 0.005656566
tailoredGlasso::sparsity(theta2.est.4!=0)
# 0.005050505
tailoredGlasso::sparsity(theta3.est.4!=0)
# 0.005050505
tailoredGlasso::sparsity(theta4.est.4!=0)
# 0.005252525

tailoredGlasso::precision(as.matrix(theta.true.4.1!=0), theta1.est.4!=0)
# 0.9642857
tailoredGlasso::precision(as.matrix(theta.true.4.2!=0), theta2.est.4!=0)
# 0.96
tailoredGlasso::precision(as.matrix(theta.true.4.3!=0), theta3.est.4!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.4.4!=0), theta4.est.4!=0)
# 1

tailoredGlasso::recall(as.matrix(theta.true.4.1!=0), theta1.est.4!=0)
# 0.2727273
tailoredGlasso::recall(as.matrix(theta.true.4.2!=0), theta2.est.4!=0)
# 0.2424242
tailoredGlasso::recall(as.matrix(theta.true.4.3!=0), theta3.est.4!=0)
# 0.2525253
tailoredGlasso::recall(as.matrix(theta.true.4.4!=0), theta4.est.4!=0)
# 0.2626263


# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.4.1 <- fastGHS(x.sf.scaled.4.1,tau_sq = 0.013,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.1 <- cov2cor(res.ecmghs.4.1$theta)
theta.est.ecmghs.4.1[which(abs(theta.est.ecmghs.4.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.1!=0)
# 0.005656566
tailoredGlasso::precision(as.matrix(theta.true.4.1!=0), theta.est.ecmghs.4.1!=0)
# 0.8571429
tailoredGlasso::recall(as.matrix(theta.true.4.1!=0), theta.est.ecmghs.4.1!=0)
# 0.2424242

# Worse than jointGHS

set.seed(22)
res.ecmghs.4.2 <- fastGHS(x.sf.scaled.4.2,tau_sq = 0.005,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.2 <- cov2cor(res.ecmghs.4.2$theta)
theta.est.ecmghs.4.2[which(abs(theta.est.ecmghs.4.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.2!=0)
# 0.005050505
tailoredGlasso::precision(as.matrix(theta.true.4.2!=0), theta.est.ecmghs.4.2!=0)
# 0.96
tailoredGlasso::recall(as.matrix(theta.true.4.2!=0), theta.est.ecmghs.4.2!=0)
# 0.2424242

# Same as jointGHS

set.seed(22)
res.ecmghs.4.3 <- fastGHS(x.sf.scaled.4.3,tau_sq = 0.008,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.3 <- cov2cor(res.ecmghs.4.3$theta)
theta.est.ecmghs.4.3[which(abs(theta.est.ecmghs.4.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.3!=0)
# 0.005050505
tailoredGlasso::precision(as.matrix(theta.true.4.3!=0), theta.est.ecmghs.4.3!=0)
# 0.96
tailoredGlasso::recall(as.matrix(theta.true.4.3!=0), theta.est.ecmghs.4.3!=0)
# 0.2424242

# Worse than jointGHS

set.seed(22)
res.ecmghs.4.4 <- fastGHS(x.sf.scaled.4.4,tau_sq = 0.0059,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.4 <- cov2cor(res.ecmghs.4.4$theta)
theta.est.ecmghs.4.4[which(abs(theta.est.ecmghs.4.4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.4!=0)
# 0.005252525
tailoredGlasso::precision(as.matrix(theta.true.4.4!=0), theta.est.ecmghs.4.4!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.4.4!=0), theta.est.ecmghs.4.4!=0)
# 0.2626263

# Same as jointGHS

# All networks have either improved estimates in jointGHS, or the same as the single-network version. 


# Compare to result from using only two networks: 

# Use jointGHS on two of the data sets

# Select tau_sq to get similar sparsity to the jointGHS results, to get comparable results
set.seed(1234)
res.joint.4.2 = jointGHS::jointGHS(list(x.sf.scaled.4.3, x.sf.scaled.4.3), tau_sq=c(0.1, 0.1),epsilon = 1e-3, fix_tau=TRUE)

theta3.est.4.2 <- cov2cor(res.joint.4.2$theta[[1]])
theta4.est.4.2 <- cov2cor(res.joint.4.2$theta[[2]])
theta3.est.4.2[which(abs(theta3.est.4.2) < 1e-5, arr.ind = T)] = 0
theta4.est.4.2[which(abs(theta4.est.4.2) < 1e-5, arr.ind = T)] = 0

tailoredGlasso::sparsity(theta3.est.4.2!=0)
# 0.005050505
tailoredGlasso::sparsity(theta4.est.4.2!=0)
# 0.005050505

tailoredGlasso::precision(as.matrix(theta.true.4.3!=0), theta3.est.4.2!=0)
# 0.96
tailoredGlasso::precision(as.matrix(theta.true.4.4!=0), theta4.est.4.2!=0)
# 0.8
tailoredGlasso::recall(as.matrix(theta.true.4.3!=0), theta3.est.4.2!=0)
# 0.2424242
tailoredGlasso::recall(as.matrix(theta.true.4.4!=0), theta4.est.4.2!=0)
# 0.2020202

# Precision and recall is worse than jointGHS. Thus, using all four networks improves the estimates. 

# EXAMPLE 5: ten datasets from related distributions, more high-dimensional ------------------------------------------------------------------

# 10% edge disagreement with the first network in the underlying edges
# Note: the 90% agreement is only with graph 1. I.e. the other networks are not as similar to each other. 

# First network: n=100, p=100
n.5.1=100
p.5.1=100
set.seed(12345)
data.sf.5.1= huge::huge.generator(n=n.5.1, d=p.5.1,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.5.1 = data.sf.5.1$theta # True adjacency matrix
theta.true.5.1 = data.sf.5.1$omega # The precision matrix
theta.true.5.1[which(theta.true.5.1<10e-5,arr.ind=T)]=0  
g.sf.5.1=graph.adjacency(data.sf.5.1$theta,mode="undirected",diag=F) # true igraph object
x.sf.5.1 = data.sf.5.1$data # Observed attributes. nxp matrix.
x.sf.scaled.5.1= scale(x.sf.5.1) # Scale columns/variables.
s.sf.scaled.5.1 = cov(x.sf.scaled.5.1) # Empirical covariance matrix
data.sf.5.1$sparsity # True sparsity: 0.02

K.5 = 10
n.vals.5 = c(n.5.1, 200, 150, 100, 250, 100, 200, 150, 100, 250)
X.list.k = list(x.sf.scaled.5.1)
theta.true.list.5 = list(theta.true.5.1)
set.seed(1234)
for (k in 2:K.5){
  graph.k = mutate.graph(data.sf.5.1, 0.1)
  X.list.k[[k]] = graph.k$data
  theta.true.list.5[[k]] = graph.k$prec.mat
}
unlist(lapply(theta.true.list.5, tailoredGlasso::sparsity))
# 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02

# Use jointGHS on the ten related data sets
set.seed(1234)
res.joint.5 = jointGHS::jointGHS(X.list.k, tau_sq=rep(1,K.5),epsilon = 1e-3, fix_tau=TRUE)
thetas.est.5 = lapply(res.joint.5$theta, cov2cor)
for(k in 1:K.5){
  thetas.est.5[[k]][which(abs(thetas.est.5[[k]]) < 1e-5, arr.ind = T)] = 0
}

# Sparsities
unlist(lapply(thetas.est.5, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
# 0.005454545 0.006464646 0.005050505 0.004848485 0.005656566 0.004444444 0.005454545 0.005454545 0.005656566 0.005252525

# Precisions
unlist(lapply(1:K.5, FUN = function(k) tailoredGlasso::precision(theta.true.list.5[[k]]!=0, thetas.est.5[[k]]!=0)))
# 0.9629630 0.7812500 0.9200000 0.9166667 0.9285714 0.9545455 0.9629630 0.9629630 0.8571429 0.9615385

# Recalls
unlist(lapply(1:K.5, FUN = function(k) tailoredGlasso::recall(theta.true.list.5[[k]]!=0, thetas.est.5[[k]]!=0)))
# 0.2626263 0.2525253 0.2323232 0.2222222 0.2626263 0.2121212 0.2626263 0.2626263 0.2424242 0.2525253

# Still feasible to compute. 

# Compare to three-network version
set.seed(1234)
res.joint.5.2 = jointGHS::jointGHS(list(X.list.k[[3]], X.list.k[[4]], X.list.k[[5]]), tau_sq=rep(10,3),epsilon = 1e-3, fix_tau=TRUE)
thetas.est.5.2 = lapply(res.joint.5.2$theta, cov2cor)
for(k in 1:3){
  thetas.est.5.2[[k]][which(abs(thetas.est.5.2[[k]]) < 1e-5, arr.ind = T)] = 0
}

# Sparsities
unlist(lapply(thetas.est.5.2, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
# 0.004444444 0.003838384 0.004848485
# Ten-network version: 0.005050505 0.004848485 0.005656566 

theta.true.list.5.2 = list(theta.true.list.5[[3]], theta.true.list.5[[4]], theta.true.list.5[[5]])
# Precisions
unlist(lapply(1:3, FUN = function(k) tailoredGlasso::precision(theta.true.list.5.2[[k]]!=0, thetas.est.5.2[[k]]!=0)))
# 0.8636364 1.0000000 0.8750000
# Ten-network version: 0.9200000 0.9166667 0.9285714 

# Recalls
unlist(lapply(1:3, FUN = function(k) tailoredGlasso::recall(theta.true.list.5.2[[k]]!=0, thetas.est.5.2[[k]]!=0)))
# 0.1919192 0.1919192 0.2121212
# Ten-network version: 0.2323232 0.2222222 0.2626263

# As we see, using more networks with jointGHS gives better results for all networks, expect wrt the precision of network 4 (as it is much sparser in the three-network version.). 




# EXAMPLE 6: look at how the sparsity changes with tau  ------------------------------------------------------------------

# 20% edge disagreement in the underlying edges

# First network: n=100, p=50, larger partial correlations (0.229)
n.6.1=100
p.6.1=50
set.seed(12345)
data.sf.6.1= huge::huge.generator(n=n.6.1, d=p.6.1,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.6.1 = data.sf.6.1$theta # True adjacency matrix
theta.true.6.1 = data.sf.6.1$omega # The precision matrix
theta.true.6.1[which(theta.true.6.1<10e-5,arr.ind=T)]=0  
g.sf.6.1=graph.adjacency(data.sf.6.1$theta,mode="undirected",diag=F) # true igraph object
x.sf.6.1 = data.sf.6.1$data # Observed attributes. nxp matrix.
x.sf.scaled.6.1= scale(x.sf.6.1) # Scale columns/variables.
s.sf.scaled.6.1 = cov(x.sf.scaled.6.1) # Empirical covariance matrix
data.sf.6.1$sparsity # True sparsity: 0.04

# Generate second data set with same sparsity
# n=200, p=50, 20% edge disagreement
n.6.2=200
p.6.2=50
set.seed(123456)
graph.6.2 = mutate.graph(data.sf.6.1, 0.2)
theta.true.6.2 = graph.6.2$prec.mat
tailoredGlasso::sparsity(theta.true.6.2!=0)
# 0.04
x.sf.scaled.6.2 = scale(mvtnorm::rmvnorm(n.6.2, sigma = solve(graph.6.2$prec.mat)))

# Use jointGHS on two slightly related data sets, with various tau values

tau.vals.6 = seq(1e-5, 10, by = 0.001)
res.joint.6 = lapply(tau.vals.6, FUN = function(t) jointGHS::jointGHS(list(x.sf.scaled.6.1, x.sf.scaled.6.2), tau_sq=c(t, t),epsilon = 1e-3, fix_tau=TRUE))
thetas.est.6.1 = lapply(res.joint.6, FUN = function(s) cov2cor(s$theta[[1]]))
thetas.est.6.2 = lapply(res.joint.6, FUN = function(s) cov2cor(s$theta[[2]]))
for(i in 1:length(tau.vals.6)){
  thetas.est.6.1[[i]][which(abs(thetas.est.6.1[[i]]) < 1e-5, arr.ind = T)] = 0
  thetas.est.6.2[[i]][which(abs(thetas.est.6.2[[i]]) < 1e-5, arr.ind = T)] = 0
}
# Sparsities
spars.6.1 = unlist(lapply(thetas.est.6.1, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
spars.6.2 = unlist(lapply(thetas.est.6.2, FUN = function(s) tailoredGlasso::sparsity(s!=0)))

data.6.1 = data.frame(tau=tau.vals.6, sparsity=spars.6.1)
data.6.2 = data.frame(tau=tau.vals.6, sparsity=spars.6.2)
p1 <- ggplot2::ggplot(data.6.1, aes(x=tau,y=sparsity))+ geom_point(color='steelblue')+ labs(title="Graph 1")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot2::ggplot(data.6.2, aes(x=tau,y=sparsity))+ geom_point(color='steelblue')+ labs(title="Graph 2")+theme(plot.title = element_text(hjust = 0.5))

pdf('examples/sparsity_vs_tau.pdf')
gridExtra::grid.arrange(p1, p2, nrow=1)
dev.off()

# As we see, the sparsity stabilized very quickly. Must look at smaller values of tau (<0.005). 

# Try grid of smaller tau values

tau.vals.6.2 = seq(1e-6, 2e-3, by = 1e-6)
res.joint.6.2 = lapply(tau.vals.6.2, FUN = function(t) jointGHS::jointGHS(list(x.sf.scaled.6.1, x.sf.scaled.6.2), tau_sq=c(t, t),epsilon = 1e-3, fix_tau=TRUE))
thetas.est.6.1.2 = lapply(res.joint.6.2, FUN = function(s) cov2cor(s$theta[[1]]))
thetas.est.6.2.2 = lapply(res.joint.6.2, FUN = function(s) cov2cor(s$theta[[2]]))
for(i in 1:length(tau.vals.6.2)){
  thetas.est.6.1.2[[i]][which(abs(thetas.est.6.1.2[[i]]) < 1e-5, arr.ind = T)] = 0
  thetas.est.6.2.2[[i]][which(abs(thetas.est.6.2.2[[i]]) < 1e-5, arr.ind = T)] = 0
}
# Sparsities
spars.6.1.2 = unlist(lapply(thetas.est.6.1.2, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
spars.6.2.2 = unlist(lapply(thetas.est.6.2.2, FUN = function(s) tailoredGlasso::sparsity(s!=0)))

# Precisions 
prec.6.1.2 = unlist(lapply(thetas.est.6.1.2, FUN = function(s) tailoredGlasso::precision(theta.true.6.1!=0, s!=0)))
prec.6.2.2 = unlist(lapply(thetas.est.6.2.2, FUN = function(s) tailoredGlasso::precision(theta.true.6.2!=0, s!=0)))

# Recalls 
rec.6.1.2 = unlist(lapply(thetas.est.6.1.2, FUN = function(s) tailoredGlasso::recall(theta.true.6.1!=0, s!=0)))
rec.6.2.2 = unlist(lapply(thetas.est.6.2.2, FUN = function(s) tailoredGlasso::recall(theta.true.6.2!=0, s!=0)))

data.6.1.2 = data.frame(tau=tau.vals.6.2, sparsity=spars.6.1.2)
data.6.2.2 = data.frame(tau=tau.vals.6.2, sparsity=spars.6.2.2)
p1 <- ggplot2::ggplot(data.6.1.2, aes(x=tau,y=sparsity))+ geom_point(color='steelblue')+ labs(title="Graph 1")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot2::ggplot(data.6.2.2, aes(x=tau,y=sparsity))+ geom_point(color='steelblue')+ labs(title="Graph 2")+theme(plot.title = element_text(hjust = 0.5))

data.prec.6.1.2 = data.frame(tau=tau.vals.6.2, precision=prec.6.1.2)
data.prec.6.2.2 = data.frame(tau=tau.vals.6.2, precision=prec.6.2.2)
p1.prec <- ggplot2::ggplot(data.prec.6.1.2, aes(x=tau,y=precision))+ geom_point(color='steelblue')+labs(title="Graph 1")+theme(plot.title = element_text(hjust = 0.5))
p2.prec <- ggplot2::ggplot(data.prec.6.2.2, aes(x=tau,y=precision))+ geom_point(color='steelblue')+labs(title="Graph 2")+theme(plot.title = element_text(hjust = 0.5))

data.rec.6.1.2 = data.frame(tau=tau.vals.6.2, recall=rec.6.1.2)
data.rec.6.2.2 = data.frame(tau=tau.vals.6.2, recall=rec.6.2.2)
p1.rec <- ggplot2::ggplot(data.rec.6.1.2, aes(x=tau,y=recall))+ geom_point(color='steelblue')+labs(title="Graph 1")+theme(plot.title = element_text(hjust = 0.5))
p2.rec <- ggplot2::ggplot(data.rec.6.2.2, aes(x=tau,y=recall))+ geom_point(color='steelblue')+labs(title="Graph 2")+theme(plot.title = element_text(hjust = 0.5))

# Plot sparsity as a function of tau

pdf('examples/sparsity_vs_tau_smallergrid.pdf')
gridExtra::grid.arrange(p1, p2, nrow=1)
dev.off()


# Plot sparsity, precision and recall as functions of tau

pdf('examples/measures_vs_tau.pdf', 10, 10)
gridExtra::grid.arrange(p1, p2, p1.prec, p2.prec, p1.rec, p2.rec,nrow=3)
dev.off()







