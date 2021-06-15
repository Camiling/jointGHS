library(fastGHS)
library(huge)
library(glasso)
library(igraph)
source('GHS.R')

## See how we can get the correct magnitude of precision matrix elements by fixing tau

# Note that now thresholding is not an issue: we either get precision matrix elements of the correct magnitude (0.1-0.2), 
# or almost zero (<1e-5) => easy to theshold and perform valiable selection. 

# Exact value of tau is not too important => a wide range can give the same sparsity level. Lambda adapts to its size. 


# EXAMPLE 1 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=50, larger partial correlations (0.229)

n.fix=100
p.fix=50
set.seed(12345)
data.sf.fix = huge::huge.generator(n=n.fix, d=p.fix,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix = data.sf.fix$theta # True adjacency matrix
theta.true.fix = data.sf.fix$omega # The precision matrix
theta.true.fix[which(theta.true.fix<10e-5,arr.ind=T)]=0  
g.sf.fix=graph.adjacency(data.sf.fix$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix = data.sf.fix$data # Observed attributes. nxp matrix.
x.sf.scaled.fix= scale(x.sf.fix) # Scale columns/variables.
s.sf.scaled.fix = cov(x.sf.scaled.fix) # Empirical covariance matrix
data.sf.fix$sparsity # True sparsity: 0.04
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2290578 0.0000000 0.0000000 0.0000000
#[2,] 0.2290578 1.0000000 0.2290578 0.0000000 0.0000000
#[3,] 0.0000000 0.2290578 1.0000000 0.2290578 0.2290578
#[4,] 0.0000000 0.0000000 0.2290578 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2290578 0.0000000 1.0000000

# Ordinary GHS
set.seed(123)
ghs.res.fix = GHS(t(x.sf.fix)%*%x.sf.fix,n.fix,burnin=100,nmc=1000)
hist(ghs.res.fix$taus.samples) # Small tau, around size 1e-4
theta.est.ghs = cov2cor(apply(ghs.res.fix$thetas.sampled, c(1,2), mean))
theta.est.off.diag.ghs <- theta.est.ghs
diag(theta.est.off.diag.ghs) <- NA
# Look at the distribution of the precision matrix elements
hist(c(theta.est.off.diag.ghs), breaks=100)
theta.est.ghs[which(abs(theta.est.ghs) < quantile(abs(theta.est.off.diag.ghs), 0.96,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.ghs), 0.96,na.rm = T)
# 0.05488769 
mean(ghs.res.fix$taus.samples)
# 0.0001879007
tailoredGlasso::sparsity(theta.est.ghs!=0)
# 0.04
tailoredGlasso::precision(theta.true.fix!=0,theta.est.ghs!=0)
# 0.7142857

theta.est.ghs[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1403512 0.0000000 0.0000000 0.0000000
#[2,] 0.1403512 1.0000000 0.2281446 0.0000000 0.0000000
#[3,] 0.0000000 0.2281446 1.0000000 0.1833029 0.2419605
#[4,] 0.0000000 0.0000000 0.1833029 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2419605 0.0000000 1.0000000

# ECMGHS

res.fix <- fastGHS(x.sf.fix,tau_sq = 0.1,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix <- cov2cor(res.fix$theta)
theta.est.fix[which(abs(theta.est.fix) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix!=0)
# 0.01632653
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 1
theta.est.fix[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1965087 0.0000000 0.0000000 0.0000000
#[2,] 0.1965087 1.0000000 0.2468123 0.0000000 0.0000000
#[3,] 0.0000000 0.2468123 1.0000000 0.2176918 0.2528399
#[4,] 0.0000000 0.0000000 0.2176918 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2528399 0.0000000 1.0000000



# Notably, we need a larger tau in ECMGHS than MCMC GHS. 



# EXAMPLE 2 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=200, p=150, larger partial correlations (0.183)

n.fix2=200
p.fix2=150
set.seed(12345)
data.sf.fix2 = huge::huge.generator(n=n.fix2, d=p.fix2,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix2 = data.sf.fix2$theta # True adjacency matrix
theta.true.fix2 = data.sf.fix2$omega # The precision matrix
theta.true.fix2[which(theta.true.fix2<10e-5,arr.ind=T)]=0  
g.sf.fix2=graph.adjacency(data.sf.fix2$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix2 = data.sf.fix2$data # Observed attributes. nxp matrix.
x.sf.scaled.fix2= scale(x.sf.fix2) # Scale columns/variables.
s.sf.scaled.fix2 = cov(x.sf.scaled.fix2) # Empirical covariance matrix
data.sf.fix2$sparsity # True sparsity: 0.013333
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix2[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1826553 0.0000000 0.0000000 0.0000000
#[2,] 0.1826553 1.0000000 0.1826553 0.0000000 0.0000000
#[3,] 0.0000000 0.1826553 1.0000000 0.1826553 0.1826553
#[4,] 0.0000000 0.0000000 0.1826553 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1826553 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here


res.fix2 <- fastGHS(x.sf.fix2,tau_sq = 0.1,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix2 <- cov2cor(res.fix2$theta)
theta.est.fix2[which(abs(theta.est.fix2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix2!=0)
# 0.007248322
tailoredGlasso::precision(as.matrix(theta.true.fix2!=0), theta.est.fix2!=0)
# 0.9012346
theta.est.fix2[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1918803 0.0000000 0.0000000 0.0000000
#[2,] 0.1918803 1.0000000 0.2173916 0.0000000 0.0000000
#[3,] 0.0000000 0.2173916 1.0000000 0.2154476 0.1886733
#[4,] 0.0000000 0.0000000 0.2154476 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1886733 0.0000000 1.0000000

# EXAMPLE 3 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=300, p=400, smaller partial correlations (0.147)

n.fix3=300
p.fix3=400
set.seed(12345)
data.sf.fix3 = huge::huge.generator(n=n.fix3, d=p.fix3,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix3 = data.sf.fix3$theta # True adjacency matrix
theta.true.fix3 = data.sf.fix3$omega # The precision matrix
theta.true.fix3[which(theta.true.fix3<10e-5,arr.ind=T)]=0  
g.sf.fix3=graph.adjacency(data.sf.fix3$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix3 = data.sf.fix3$data # Observed attributes. nxp matrix.
x.sf.scaled.fix3= scale(x.sf.fix3) # Scale columns/variables.
s.sf.scaled.fix3 = cov(x.sf.scaled.fix3) # Empirical covariance matrix
data.sf.fix3$sparsity # True sparsity: 0.005
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix3[1:5,1:5])
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1530514 0.0000000 0.0000000 0.0000000
#[2,] 0.1530514 1.0000000 0.1530514 0.0000000 0.0000000
#[3,] 0.0000000 0.1530514 1.0000000 0.1530514 0.1530514
#[4,] 0.0000000 0.0000000 0.1530514 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1530514 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here

res.fix3 <- fastGHS(x.sf.fix3,tau_sq = 0.01,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix3 <- cov2cor(res.fix3$theta)
theta.est.fix3[which(abs(theta.est.fix3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix3!=0)
# 0.006704261
tailoredGlasso::precision(as.matrix(theta.true.fix3!=0), theta.est.fix3!=0)
# 0.4168224
theta.est.fix3[1:5,1:5]
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1901565 0.0000000 0.0000000 0.0000000
#[2,] 0.1901565 1.0000000 0.1400074 0.0000000 0.0000000
#[3,] 0.0000000 0.1400074 1.0000000 0.1752922 0.1859689
#[4,] 0.0000000 0.0000000 0.1752922 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1859689 0.0000000 1.0000000

# EXAMPLE 4 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=200, smaller partial correlations (0.18)

n.fix4=100
p.fix4=200
set.seed(12345)
data.sf.fix4 = huge::huge.generator(n=n.fix4, d=p.fix4,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix4 = data.sf.fix4$theta # True adjacency matrix
theta.true.fix4 = data.sf.fix4$omega # The precision matrix
theta.true.fix4[which(theta.true.fix4<10e-5,arr.ind=T)]=0  
g.sf.fix4=graph.adjacency(data.sf.fix4$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix4 = data.sf.fix4$data # Observed attributes. nxp matrix.
x.sf.scaled.fix4= scale(x.sf.fix4) # Scale columns/variables.
s.sf.scaled.fix4 = cov(x.sf.scaled.fix4) # Empirical covariance matrix
data.sf.fix4$sparsity # True sparsity: 0.01
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix4[1:5,1:5])
#           [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1807393 0.0000000 0.0000000 0.0000000
#[2,] 0.1807393 1.0000000 0.1807393 0.0000000 0.0000000
#[3,] 0.0000000 0.1807393 1.0000000 0.1807393 0.1807393
#[4,] 0.0000000 0.0000000 0.1807393 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1807393 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here

res.fix4 <- fastGHS(x.sf.fix4,tau_sq = 0.025,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix4 <- cov2cor(res.fix4$theta)
theta.est.fix4[which(abs(theta.est.fix4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix4!=0)
# 0.01005025
tailoredGlasso::precision(as.matrix(theta.true.fix4!=0), theta.est.fix4!=0)
# 0.355
theta.est.fix4[1:5,1:5]
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000 0.000000 0.0000000
#[2,]    0 1.0000000 0.2065228 0.000000 0.0000000
#[3,]    0 0.2065228 1.0000000 0.186527 0.2750568
#[4,]    0 0.0000000 0.1865270 1.000000 0.0000000
#[5,]    0 0.0000000 0.2750568 0.000000 1.0000000

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate

gg = huge(x.sf.fix4,method='glasso', lambda=0.2879)
gg$sparsity
# 0.01005025
tailoredGlasso::precision(as.matrix(theta.true.fix4!=0), gg$icov[[1]]!=0)
# 0.35

# Better precision for ECM GHS than the graphical lasso!



