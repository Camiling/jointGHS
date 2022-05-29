library(jointGHS)
library(huge)
library(glasso)
library(igraph)
library(fastGHS)
library(ggplot2)
library(gridExtra)
library(foreach)
library(tailoredGlasso) # Install from Github/Camiling/tailoredGlasso
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
res.joint.1 = jointGHS::jointGHS(list(x.sf.scaled.1, x2.sf.1.scaled), epsilon = 1e-5, AIC_selection = T, AIC_eps = 0.1)

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
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000 0.000000 0.0000000
#[2,]    0 1.0000000 0.3312751 0.000000 0.0000000
#[3,]    0 0.3312751 1.0000000 0.246183 0.2389791
#[4,]    0 0.0000000 0.2461830 1.000000 0.0000000
#[5,]    0 0.0000000 0.2389791 0.000000 1.0000000


res.joint.1$tau_sq
# 6.201 5.201

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.1 <- fastGHS(x.sf.scaled.1,tau_sq = 0.02,epsilon = 1e-5, fix_tau=TRUE)
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
res.joint.2 = jointGHS::jointGHS(list(x.sf.scaled.2.1, x.sf.scaled.2.2), AIC_selection = T, epsilon = 1e-3, AIC_eps = 0.01)

theta1.est.2 <- cov2cor(res.joint.2$theta[[1]])
theta2.est.2 <- cov2cor(res.joint.2$theta[[2]])
theta1.est.2[which(abs(theta1.est.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta1.est.2!=0)
# 0.008163265
theta2.est.2[which(abs(theta2.est.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta2.est.2!=0)
# 0.005714286
tailoredGlasso::precision(as.matrix(theta.true.2.1!=0), theta1.est.2!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.2.2!=0), theta2.est.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.2.1!=0), theta1.est.2!=0)
# 0.2040816
tailoredGlasso::recall(as.matrix(theta.true.2.2!=0), theta2.est.2!=0)
# 0.1428571

res.joint.2$tau_sq
# 9.601 7.201

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.2.1 <- fastGHS(x.sf.scaled.2.1,tau_sq = 0.0082,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.2.1 <- cov2cor(res.ecmghs.2.1$theta)
theta.est.ecmghs.2.1[which(abs(theta.est.ecmghs.2.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.2.1!=0)
# 0.008163265
tailoredGlasso::precision(as.matrix(theta.true.2.1!=0), theta.est.ecmghs.2.1!=0)
# 0.9
tailoredGlasso::recall(as.matrix(theta.true.2.1!=0), theta.est.ecmghs.2.1!=0)
# 0.1836735

set.seed(22)
res.ecmghs.2.2 <- fastGHS(x.sf.scaled.2.2,tau_sq = 0.0025,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.2.2 <- cov2cor(res.ecmghs.2.2$theta)
theta.est.ecmghs.2.2[which(abs(theta.est.ecmghs.2.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.2.2!=0)
# 0.005714286
tailoredGlasso::precision(as.matrix(theta.true.2.2!=0), theta.est.ecmghs.2.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.2.2!=0), theta.est.ecmghs.2.2!=0)
# 0.1428571

# For competely unrelated networks, jointGHS gives approximately the same results as the ordinary ECM GHS 
# (same for the network with large n, better for the one with the smaller n)

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
set.seed(1234567)
graph.3.2 = mutate.graph.old(data.sf.3.1, 0.2)
theta.true.3.2 = graph.3.2$prec.mat
tailoredGlasso::sparsity(theta.true.3.2!=0)
# 0.04
x.sf.scaled.3.2 = scale(mvtnorm::rmvnorm(n.3.2, sigma = solve(graph.3.2$prec.mat)))

# Use jointGHS on two slightly related data sets
res.joint.3 = jointGHS::jointGHS(list(x.sf.scaled.3.1, x.sf.scaled.3.2), AIC_selection = T, epsilon = 1e-3, AIC_eps = 0.1)

theta1.est.3 <- cov2cor(res.joint.3$theta[[1]])
theta2.est.3 <- cov2cor(res.joint.3$theta[[2]])
theta1.est.3[which(abs(theta1.est.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta1.est.3!=0)
# 0.01387755
theta2.est.3[which(abs(theta2.est.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta2.est.3!=0)
# 0.01387755
tailoredGlasso::precision(as.matrix(theta.true.3.1!=0), theta1.est.3!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.3.2!=0), theta2.est.3!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.3.1!=0), theta1.est.3!=0)
# 0.3469388
tailoredGlasso::recall(as.matrix(theta.true.3.2!=0), theta2.est.3!=0)
# 0.3469388

res.joint.3$tau_sq
# 6.001 5.001

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.3.1 <- fastGHS(x.sf.scaled.3.1,tau_sq = 0.025,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.3.1 <- cov2cor(res.ecmghs.3.1$theta)
theta.est.ecmghs.3.1[which(abs(theta.est.ecmghs.3.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.3.1!=0)
# 0.01387755
tailoredGlasso::precision(as.matrix(theta.true.3.1!=0), theta.est.ecmghs.3.1!=0)
# 0.8823529
tailoredGlasso::recall(as.matrix(theta.true.3.1!=0), theta.est.ecmghs.3.1!=0)
# 0.3061224

set.seed(22)
res.ecmghs.3.2 <- fastGHS(x.sf.scaled.3.2,tau_sq = 0.0032,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.3.2 <- cov2cor(res.ecmghs.3.2$theta)
theta.est.ecmghs.3.2[which(abs(theta.est.ecmghs.3.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.3.2!=0)
# 0.01387755
tailoredGlasso::precision(as.matrix(theta.true.3.2!=0), theta.est.ecmghs.3.2!=0)
# 0.9411765
tailoredGlasso::recall(as.matrix(theta.true.3.2!=0), theta.est.ecmghs.3.2!=0)
# 0.3265306

# For related networks, jointGHS gives more accurate results than the single network ECM GHS 


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
graph.4.2 = mutate.graph.old(data.sf.4.1, 0.2)
theta.true.4.2 = graph.4.2$prec.mat
tailoredGlasso::sparsity(theta.true.4.2!=0)
# 0.02
x.sf.scaled.4.2 = scale(mvtnorm::rmvnorm(n.4.2, sigma = solve(graph.4.2$prec.mat)))

# Generate third data set with same sparsity
# n=150, p=100, 20% edge disagreement
n.4.3=150
p.4.3=100
set.seed(12345677)
graph.4.3 = mutate.graph.old(data.sf.4.1, 0.2)
theta.true.4.3 = graph.4.3$prec.mat
tailoredGlasso::sparsity(theta.true.4.3!=0)
# 0.02
x.sf.scaled.4.3 = scale(mvtnorm::rmvnorm(n.4.3, sigma = solve(graph.4.3$prec.mat)))

# Generate fourth data set with same sparsity
# n=250, p=100, 20% edge disagreement
n.4.4=250
p.4.4=100
set.seed(12345688)
graph.4.4 = mutate.graph.old(data.sf.4.1, 0.2)
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
res.joint.4 = jointGHS::jointGHS(list(x.sf.scaled.4.1, x.sf.scaled.4.2, x.sf.scaled.4.3, x.sf.scaled.4.4), AIC_selection = T, epsilon = 1e-3, AIC_eps = 0.1)

theta1.est.4 <- cov2cor(res.joint.4$theta[[1]])
theta2.est.4 <- cov2cor(res.joint.4$theta[[2]])
theta3.est.4 <- cov2cor(res.joint.4$theta[[3]])
theta4.est.4 <- cov2cor(res.joint.4$theta[[4]])
theta1.est.4[which(abs(theta1.est.4) < 1e-5, arr.ind = T)] = 0
theta2.est.4[which(abs(theta2.est.4) < 1e-5, arr.ind = T)] = 0
theta3.est.4[which(abs(theta3.est.4) < 1e-5, arr.ind = T)] = 0
theta4.est.4[which(abs(theta4.est.4) < 1e-5, arr.ind = T)] = 0

tailoredGlasso::sparsity(theta1.est.4!=0)
# 0.005050505
tailoredGlasso::sparsity(theta2.est.4!=0)
# 0.004646465
tailoredGlasso::sparsity(theta3.est.4!=0)
# 0.004444444
tailoredGlasso::sparsity(theta4.est.4!=0)
# 0.004848485

tailoredGlasso::precision(as.matrix(theta.true.4.1!=0), theta1.est.4!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.4.2!=0), theta2.est.4!=0)
# 0.9565217
tailoredGlasso::precision(as.matrix(theta.true.4.3!=0), theta3.est.4!=0)
# 1
tailoredGlasso::precision(as.matrix(theta.true.4.4!=0), theta4.est.4!=0)
# 1

tailoredGlasso::recall(as.matrix(theta.true.4.1!=0), theta1.est.4!=0)
# 0.2525253
tailoredGlasso::recall(as.matrix(theta.true.4.2!=0), theta2.est.4!=0)
# 0.2222222
tailoredGlasso::recall(as.matrix(theta.true.4.3!=0), theta3.est.4!=0)
# 0.2222222
tailoredGlasso::recall(as.matrix(theta.true.4.4!=0), theta4.est.4!=0)
# 0.2424242

res.joint.4$tau_sq
# 12.401  5.201  7.601  3.401

# Single-network ECM for GHS on each network separately, forced to the same sparsity to allow for direct comparison

set.seed(22)
res.ecmghs.4.1 <- fastGHS(x.sf.scaled.4.1,tau_sq = 0.013,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.1 <- cov2cor(res.ecmghs.4.1$theta)
theta.est.ecmghs.4.1[which(abs(theta.est.ecmghs.4.1) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.1!=0)
# 0.004646465
tailoredGlasso::precision(as.matrix(theta.true.4.1!=0), theta.est.ecmghs.4.1!=0)
# 0.9565217
tailoredGlasso::recall(as.matrix(theta.true.4.1!=0), theta.est.ecmghs.4.1!=0)
# 0.2222222

# Worse than jointGHS

set.seed(22)
res.ecmghs.4.2 <- fastGHS(x.sf.scaled.4.2,tau_sq = 0.003495,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.2 <- cov2cor(res.ecmghs.4.2$theta)
theta.est.ecmghs.4.2[which(abs(theta.est.ecmghs.4.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.2!=0)
# 0.003838384
tailoredGlasso::precision(as.matrix(theta.true.4.2!=0), theta.est.ecmghs.4.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.4.2!=0), theta.est.ecmghs.4.2!=0)
# 0.1919192

# Similar to jointGHS

set.seed(22)
res.ecmghs.4.3 <- fastGHS(x.sf.scaled.4.3,tau_sq = 0.006,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.3 <- cov2cor(res.ecmghs.4.3$theta)
theta.est.ecmghs.4.3[which(abs(theta.est.ecmghs.4.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.3!=0)
# 0.004040404
tailoredGlasso::precision(as.matrix(theta.true.4.3!=0), theta.est.ecmghs.4.3!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.4.3!=0), theta.est.ecmghs.4.3!=0)
# 0.2020202

# Worse than jointGHS

set.seed(22)
res.ecmghs.4.4 <- fastGHS(x.sf.scaled.4.4,tau_sq = 0.0025,epsilon = 1e-3, fix_tau=TRUE)
theta.est.ecmghs.4.4 <- cov2cor(res.ecmghs.4.4$theta)
theta.est.ecmghs.4.4[which(abs(theta.est.ecmghs.4.4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.ecmghs.4.4!=0)
# 0.003636364
tailoredGlasso::precision(as.matrix(theta.true.4.4!=0), theta.est.ecmghs.4.4!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.4.4!=0), theta.est.ecmghs.4.4!=0)
# 0.1818182

# Worse than jointGHS

# Better with jointGHS, biggest improvement for the network with the fewest observations


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
  graph.k = mutate.graph.old(data.sf.5.1, 0.1)
  X.list.k[[k]] = graph.k$data
  theta.true.list.5[[k]] = graph.k$prec.mat
}
unlist(lapply(theta.true.list.5, tailoredGlasso::sparsity))
# 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02

# Use jointGHS on the ten related data sets
set.seed(1234)
res.joint.5 = jointGHS::jointGHS(X.list.k, AIC_selection = T,epsilon = 1e-3, AIC_eps = 0.1)
thetas.est.5 = lapply(res.joint.5$theta, cov2cor)
for(k in 1:K.5){
  thetas.est.5[[k]][which(abs(thetas.est.5[[k]]) < 1e-5, arr.ind = T)] = 0
}

# Sparsities
unlist(lapply(thetas.est.5, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
# 0.004848485 0.006464646 0.005050505 0.004848485 0.005454545 0.004646465 0.005454545 0.005454545 0.005454545 0.005050505

# Precisions
unlist(lapply(1:K.5, FUN = function(k) tailoredGlasso::precision(theta.true.list.5[[k]]!=0, thetas.est.5[[k]]!=0)))
# 1.0000000 0.7812500 0.9200000 0.9166667 0.9259259 0.9565217 0.9629630 0.9629630 0.8518519 1.0000000

# Recalls
unlist(lapply(1:K.5, FUN = function(k) tailoredGlasso::recall(theta.true.list.5[[k]]!=0, thetas.est.5[[k]]!=0)))
# 0.2424242 0.2525253 0.2323232 0.2222222 0.2525253 0.2222222 0.2626263 0.2626263 0.2323232 0.2525253

# Computes in a few minutes

# Compare to three-network version
set.seed(1234)
res.joint.5.2 = jointGHS::jointGHS(list(X.list.k[[3]], X.list.k[[4]], X.list.k[[5]]), AIC_selection = T,epsilon = 1e-3, AIC_eps = 0.1)
thetas.est.5.2 = lapply(res.joint.5.2$theta, cov2cor)
for(k in 1:3){
  thetas.est.5.2[[k]][which(abs(thetas.est.5.2[[k]]) < 1e-5, arr.ind = T)] = 0
}

# Sparsities
unlist(lapply(thetas.est.5.2, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
# 0.004444444 0.003838384 0.004646465
# Ten-network version: 0.005050505 0.004848485 0.005454545

theta.true.list.5.2 = list(theta.true.list.5[[3]], theta.true.list.5[[4]], theta.true.list.5[[5]])
# Precisions
unlist(lapply(1:3, FUN = function(k) tailoredGlasso::precision(theta.true.list.5.2[[k]]!=0, thetas.est.5.2[[k]]!=0)))
# 0.9090909 1.0000000 0.8695652
# Ten-network version: 0.9200000 0.9166667 0.9259259

# Recalls
unlist(lapply(1:3, FUN = function(k) tailoredGlasso::recall(theta.true.list.5.2[[k]]!=0, thetas.est.5.2[[k]]!=0)))
# 0.2020202 0.1919192 0.2020202
# Ten-network version: 0.2323232 0.2222222 0.2525253

# Similar precision, better recall with more networks






