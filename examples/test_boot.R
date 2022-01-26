library(jointGHS)
library(huge)
library(glasso)
library(igraph)
library(fastGHS)
library(ggplot2)
library(gridExtra)
library(foreach)
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
set.seed(123)
res.joint.1 = jointGHS::jointGHS(list(x.sf.scaled.1, x2.sf.1.scaled), epsilon = 1e-5, AIC_selection = T, AIC_eps = 0.1, B=1000, boot_check=TRUE)

# Edge 1--2 in network 1 - a true edge, but weaker signal (not detected in joint model)
lambdas.test.k1 = unlist(lapply(res.joint.1$Lambda_sq_boot, FUN= function(s) s[[1]][1,2]))
png("examples/Figures/Lambda_sq_boot_trueweak.png", width=1200, height =1200)
hist(c(lambdas.test.k1), breaks=100, main='Lambda_sq found from bootstrap samples, for true edge with weak signal', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][1,2],2), y=c(0,300), col='red')
lines(x=rep(quantile(lambdas.test.k1, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k1, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][1,2]
quantile(lambdas.test.k1, c(0.025, 0.975))

# Edge 3--7 in network 1 - a true edge with strong signal (detected in joint model)
lambdas.test.k2 = unlist(lapply(res.joint.1$Lambda_sq_boot, FUN= function(s) s[[1]][3,7]))
png("examples/Figures/Lambda_sq_boot_truestrong.png", width=1200, height =1200)
hist(c(lambdas.test.k2), breaks=100, main='Lambda_sq found from bootstrap samples, for true edge with strong signal', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][3,7],2), y=c(0,100), col='red')
lines(x=rep(quantile(lambdas.test.k2, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k2, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][3,7]
quantile(lambdas.test.k2, c(0.025, 0.975))

# Edge 3--7 in network 1 - a true edge with strong signal (zoomed in)
lambdas.test.k2.2 = unlist(lapply(res.joint.1$Lambda_sq_boot, FUN= function(s) s[[1]][3,7]))
lambdas.test.k2.2= lambdas.test.k2[which(lambdas.test.k2.2<0.025,arr.ind=T)]
png("examples/Figures/Lambda_sq_boot_truestrong_zoomed.png", width=1200, height =1200)
hist(c(lambdas.test.k2.2), breaks=100, main='Lambda_sq found from bootstrap samples, for true edge with strong signal', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][3,7],2), y=c(0,100), col='red')
lines(x=rep(quantile(lambdas.test.k2, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k2, 0.975),2), y=c(0,100), col='blue')
dev.off()


# Edge 1--4 in network 1 - not an edge
lambdas.test.k3 = unlist(lapply(res.joint.1$Lambda_sq_boot, FUN= function(s) s[[1]][1,4]))
png("examples/Figures/Lambda_sq_boot_false.png", width=1200, height =1200)
hist(c(lambdas.test.k3), breaks=100, main='Lambda_sq found from bootstrap samples, for nonexistent edge', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][1,4],2), y=c(0,300), col='red')
lines(x=rep(quantile(lambdas.test.k3, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k3, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][1,4]
quantile(lambdas.test.k3, c(0.025, 0.975))

# Edge 2--3 in network 1 - very strong edge (detected in joint model)
lambdas.test.k4 = unlist(lapply(res.joint.1$Lambda_sq_boot, FUN= function(s) s[[1]][2,3]))
png("examples/Figures/Lambda_sq_trueverystrong.png", width=1200, height =1200)
hist(c(lambdas.test.k4), breaks=100, main='Lambda_sq found from bootstrap samples, for true edge with very strong signal', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][2,3],2), y=c(0,100), col='red')
lines(x=rep(quantile(lambdas.test.k4, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k4, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][2,3]
quantile(lambdas.test.k4, c(0.025, 0.975))


# Edge 1--13 in network 1 - not an edge
lambdas.test.k5 = unlist(lapply(res.joint.1$Lambda_sq_boot, FUN= function(s) s[[1]][1,13]))
png("examples/Figures/Lambda_sq_boot_false_2.png", width=1200, height =1200)
hist(c(lambdas.test.k5), breaks=100, main='Lambda_sq found from bootstrap samples, for nonexistent edge', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][1,13],2), y=c(0,300), col='red')
lines(x=rep(quantile(lambdas.test.k5, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k5, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][1,13]
quantile(lambdas.test.k5, c(0.025, 0.975))

## What if we discard the negative values for lambda_ijs with clear distributions for non-zero values?

# Must then have a rule for when we consider them to be non-zero, e.g. at least 50% of bootstrap values are non-zero
# Note that lambda_ij's tend to be either larger than 1e-5 or smaller than 1e-26 (smaller than machine eps)

# Edge 1--2
mean(lambdas.test.k1>1e-5)
# 0.25
# Not enough. 

# Edge 3--7
mean(lambdas.test.k2>1e-5)
# 0.72
# Enough. 
lambdas.test.k2.2 = lambdas.test.k2[which(lambdas.test.k2>1e-5)]
png("examples/Figures/Lambda_sq_boot_truestrong_onlypositive.png", width=1200, height =1200)
hist(c(lambdas.test.k2.2), breaks=100, main='Lambda_sq found from bootstrap samples, for true edge with strong signal', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][3,7],2), y=c(0,100), col='red')
lines(x=rep(quantile(lambdas.test.k2.2, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k2.2, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][3,7]
quantile(lambdas.test.k2.2, c(0.025, 0.975))

# Edge 2--3
mean(lambdas.test.k4>1e-5)
# 0.9725
# Edge 2--3 in network 1 - very strong edge (detected in joint model)
lambdas.test.k4.2 = lambdas.test.k4[which(lambdas.test.k4>1e-5)]
png("examples/Figures/Lambda_sq_trueverystrong_onlypositive.png", width=1200, height =1200)
hist(c(lambdas.test.k4.2), breaks=100, main='Lambda_sq found from bootstrap samples, for true edge with very strong signal', xlab='Lambda_sq')
lines(x=rep(res.joint.1$Lambda_sq[[1]][2,3],2), y=c(0,100), col='red')
lines(x=rep(quantile(lambdas.test.k4.2, 0.025),2), y=c(0,100), col='blue')
lines(x=rep(quantile(lambdas.test.k4.2, 0.975),2), y=c(0,100), col='blue')
dev.off()
res.joint.1$Lambda_sq[[1]][2,3]
quantile(lambdas.test.k4.2, c(0.025, 0.975))

# Test plotting and printing functionality 

# Plot BB results for all inferred edges in network 1
plot(res.joint.1, k=1)

# Plot BB results for specific potential edge
plot(res.joint.1, k=1, edges=c(1,13))

# Plot BB results for specific potential edges in network k
plot(res.joint.1, k=2, edges=matrix(c(1,13, 2,5, 6,7, 8,9),ncol=2))

# Plot all networks
plot(res.joint.1, plot_boot = F)

# Print BB results for network 1
print(res.joint.1, k=1)

# Return as data frame
res.df = print(res.joint.1, return_df=T)

# Print for ALL potential edges
print(res.joint.1, k=1, edges=t(combn(1:p.1, 2)))
# Too many edges, gives a messy printout... Better to stay with inferred edges only

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
set.seed(123)
res.joint.2 = jointGHS::jointGHS(list(x.sf.scaled.2.1, x.sf.scaled.2.2), AIC_selection = T, epsilon = 1e-3, AIC_eps = 0.01, B=1000, boot_check=TRUE)

# Posterior checks
plot(res.joint.2)

print(res.joint.2)


# EXAMPLE 3: four data sets - first three related and last unrelated ------------------------------------------------------------------


# First network: n=100, p=50
n.3=100
p.3=50
set.seed(12345)
data.sf.3= huge::huge.generator(n=n.3, d=p.3,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.3 = data.sf.3$theta # True adjacency matrix
theta.true.3 = data.sf.3$omega # The precision matrix
theta.true.3[which(theta.true.3<10e-5,arr.ind=T)]=0  
g.sf.3=graph.adjacency(data.sf.3$theta,mode="undirected",diag=F) # true igraph object
x.sf.3 = data.sf.3$data # Observed attributes. nxp matrix.
x.sf.scaled.3= scale(x.sf.3) # Scale columns/variables.
s.sf.scaled.3 = cov(x.sf.scaled.3) # Empirical covariance matrix
data.sf.3$sparsity # True sparsity: 0.04

# Generate second data set with same sparsity, 10% edge disagreement
n.3.2=200
p.3.2=50
set.seed(123456)
graph.3.2 = mutate.graph(data.sf.3, 0.1)
theta.true.3.2 = graph.3.2$prec.mat
tailoredGlasso::sparsity(theta.true.3.2!=0)
# 0.04
x.sf.scaled.3.2 = scale(mvtnorm::rmvnorm(n.3.2, sigma = solve(graph.3.2$prec.mat)))

# Generate third data set with same sparsity, 10% edge disagreement
n.3.3=150
p.3.3=50
set.seed(12345677)
graph.3.3 = mutate.graph(data.sf.3, 0.1)
theta.true.3.3 = graph.3.3$prec.mat
tailoredGlasso::sparsity(theta.true.3.3!=0)
# 0.04
x.sf.scaled.3.3 = scale(mvtnorm::rmvnorm(n.3.3, sigma = solve(graph.3.3$prec.mat)))

# Generate fourth data set with same sparsity, unrelated
n.3.4=150
p.3.4=50
set.seed(12345688)
data.sf.3.4 = huge::huge.generator(n=n.3.4, d=p.3.4,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.3.4 = data.sf.3.4$theta # True adjacency matrix
theta.true.3.4 = data.sf.3.4$omega # The precision matrix
theta.true.3.4[which(theta.true.3.4<10e-5,arr.ind=T)]=0  
x.sf.3.4 = data.sf.3.4$data # Observed attributes. nxp matrix.
x.sf.scaled.3.4= scale(x.sf.3.4) # Scale columns/variables.
data.sf.3$sparsity # True sparsity: 0.04

thetas.true.list.3 = list(theta.true.3, theta.true.3.2, theta.true.3.3, theta.true.3.4)

# Use jointGHS on the four related data sets
set.seed(1234)
res.joint.3 = jointGHS::jointGHS(list(x.sf.scaled.3, x.sf.scaled.3.2, x.sf.scaled.3.3, x.sf.scaled.3.4), AIC_selection = T, epsilon = 1e-3, AIC_eps = 0.1,B=400, boot_check=TRUE)

# Posterior checks
plot(res.joint.3)

print(res.joint.3)


