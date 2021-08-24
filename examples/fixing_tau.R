library(fastGHS)
library(huge)
library(glasso)
library(igraph)
source('GHS.R')

## See how we can get the correct magnitude of precision matrix elements by fixing tau

# We also test how initialisations and choice of tau affect the resulting estimate. 

# Note that now thresholding is not an issue: we either get precision matrix elements of the correct magnitude (0.1-0.2), 
# or almost zero (<1e-5) => easy to threshold and perform variable selection. MCMC GHS does not have the same property. 

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


# ECMGHS
system.time(
res.fix <- fastGHS(x.sf.fix,tau_sq = 0.1,epsilon = 1e-3, fix_tau=TRUE)
)
theta.est.fix <- cov2cor(res.fix$theta)
theta.est.fix[which(abs(theta.est.fix) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix!=0)
# 0.01632653
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 0.4081633
theta.est.fix[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1965087 0.0000000 0.0000000 0.0000000
#[2,] 0.1965087 1.0000000 0.2468123 0.0000000 0.0000000
#[3,] 0.0000000 0.2468123 1.0000000 0.2176918 0.2528399
#[4,] 0.0000000 0.0000000 0.2176918 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2528399 0.0000000 1.0000000

# Compare to ordinary GHS, forced to same sparsity to allow for direct comparison
set.seed(123)
system.time(
ghs.res.fix <- GHS(t(x.sf.fix)%*%x.sf.fix,n.fix,burnin=100,nmc=1000)
)
hist(ghs.res.fix$taus.samples) # Small tau, around size 1e-4
theta.est.ghs = cov2cor(apply(ghs.res.fix$thetas.sampled, c(1,2), mean))
theta.est.off.diag.ghs <- theta.est.ghs
diag(theta.est.off.diag.ghs) <- NA
# Look at the distribution of the precision matrix elements: no clear distinguishment between 
hist(c(theta.est.off.diag.ghs), breaks=100)
theta.est.ghs[which(abs(theta.est.ghs) < quantile(abs(theta.est.off.diag.ghs),1-tailoredGlasso::sparsity(theta.est.fix!=0),na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.ghs),1-tailoredGlasso::sparsity(theta.est.fix!=0),na.rm = T)
# 98.36735% 
# 0.1721756 
mean(ghs.res.fix$taus.samples)
# 0.0001879007
tailoredGlasso::sparsity(theta.est.ghs!=0)
#  0.01632653
tailoredGlasso::precision(theta.true.fix!=0,theta.est.ghs!=0)
# 1
tailoredGlasso::recall(theta.true.fix!=0,theta.est.ghs!=0)
# 0.4081633
theta.est.ghs[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000 0.0000000 0.0000000
#[2,]    0 1.0000000 0.2281446 0.0000000 0.0000000
#[3,]    0 0.2281446 1.0000000 0.1833029 0.2419605
#[4,]    0 0.0000000 0.1833029 1.0000000 0.0000000
#[5,]    0 0.0000000 0.2419605 0.0000000 1.0000000

# We get the same precision and recall with the methods. 

# We also see that ECMGHS is almost 100 times faster than the MCMC approach. 

# Notably, we need a larger tau in ECMGHS than MCMC GHS. 

# Look at their edge agreement
tailoredGlasso::confusion.matrix(theta.est.fix!=0,theta.est.ghs!=0)
#     [,1] [,2]
#[1,]   16    4
#[2,]    4 1201

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix,method='glasso', lambda=0.465)
gg$sparsity
# 0.01632653
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), gg$icov[[1]]!=0)
# 0.8
tailoredGlasso::recall(as.matrix(theta.true.fix!=0), gg$icov[[1]]!=0)
# 0.3265306

# Better precision and recall for ECM GHS than the graphical lasso!



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
tailoredGlasso::recall(as.matrix(theta.true.fix2!=0), theta.est.fix2!=0)
# 0.4899329
theta.est.fix2[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1918803 0.0000000 0.0000000 0.0000000
#[2,] 0.1918803 1.0000000 0.2173916 0.0000000 0.0000000
#[3,] 0.0000000 0.2173916 1.0000000 0.2154476 0.1886733
#[4,] 0.0000000 0.0000000 0.2154476 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1886733 0.0000000 1.0000000

# Notably, we need a larger tau in ECMGHS than MCMC GHS. 

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix2,method='glasso', lambda=0.2745)
gg$sparsity
# 0.007248322
tailoredGlasso::precision(as.matrix(theta.true.fix2!=0), gg$icov[[1]]!=0)
# 0.6419753
tailoredGlasso::recall(as.matrix(theta.true.fix2!=0), gg$icov[[1]]!=0)
# 0.3489933

# Better precision and recall for ECM GHS than the graphical lasso!


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
tailoredGlasso::recall(as.matrix(theta.true.fix3!=0), theta.est.fix3!=0)
# 0.5588972
theta.est.fix3[1:5,1:5]
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1901565 0.0000000 0.0000000 0.0000000
#[2,] 0.1901565 1.0000000 0.1400074 0.0000000 0.0000000
#[3,] 0.0000000 0.1400074 1.0000000 0.1752922 0.1859689
#[4,] 0.0000000 0.0000000 0.1752922 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1859689 0.0000000 1.0000000

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix3,method='glasso', lambda=0.1754)
gg$sparsity
# 0.006704261
tailoredGlasso::precision(as.matrix(theta.true.fix3!=0), gg$icov[[1]]!=0)
# 0.3850467
tailoredGlasso::recall(as.matrix(theta.true.fix3!=0), gg$icov[[1]]!=0)
# 0.5162907

# Better precision and recall for ECM GHS than the graphical lasso!

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
tailoredGlasso::recall(as.matrix(theta.true.fix4!=0), theta.est.fix4!=0)
# 0.3567839
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
tailoredGlasso::recall(as.matrix(theta.true.fix4!=0), gg$icov[[1]]!=0)
# 0.3517588

# Better precision and recall for ECM GHS than the graphical lasso!

# EXAMPLE 5 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=50, p=100, larger partial correlations (0.192)

# Test sensitivity to choice of tau (within same magnitude), with weak signal, much noise

n.fix5=50
p.fix5=100
set.seed(12345)
data.sf.fix5 = huge::huge.generator(n=n.fix5, d=p.fix5,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix5 = data.sf.fix5$theta # True adjacency matrix
theta.true.fix5 = data.sf.fix5$omega # The precision matrix
theta.true.fix5[which(theta.true.fix5<10e-5,arr.ind=T)]=0  
g.sf.fix5=graph.adjacency(data.sf.fix5$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix5 = data.sf.fix5$data # Observed attributes. nxp matrix.
x.sf.scaled.fix5= scale(x.sf.fix5) # Scale columns/variables.
s.sf.scaled.fix5 = cov(x.sf.scaled.fix5) # Empirical covariance matrix
data.sf.fix5$sparsity # True sparsity: 0.02
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix5[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1918105 0.0000000 0.0000000 0.0000000
#[2,] 0.1918105 1.0000000 0.1918105 0.0000000 0.0000000
#[3,] 0.0000000 0.1918105 1.0000000 0.1918105 0.1918105
#[4,] 0.0000000 0.0000000 0.1918105 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1918105 0.0000000 1.0000000

# First tau

res.fix5 <- fastGHS(x.sf.fix5,tau_sq = 0.1,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix5 <- cov2cor(res.fix5$theta)
theta.est.fix5[which(abs(theta.est.fix5) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5!=0)
# 0.0240404
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5!=0)
# 0.2521008
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5!=0)
# 0.3030303
theta.est.fix5[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000  0.0000000  0.0000000
#[2,]    0 1.0000000 0.2624593  0.0000000  0.0000000
#[3,]    0 0.2624593 1.0000000  0.2510143  0.0000000
#[4,]    0 0.0000000 0.2510143  1.0000000 -0.3653482
#[5,]    0 0.0000000 0.0000000 -0.3653482  1.0000000

# Second tau

res.fix5.2 <- fastGHS(x.sf.fix5,tau_sq = 0.075,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix5.2 <- cov2cor(res.fix5.2$theta)
theta.est.fix5.2[which(abs(theta.est.fix5.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.2!=0)
# 0.02080808
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.2!=0)
# 0.2621359
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.2!=0)
# 0.2727273
theta.est.fix5.2[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000  0.0000000  0.000000
#[2,]    0 1.0000000 0.2757929  0.0000000  0.000000
#[3,]    0 0.2757929 1.0000000  0.2503962  0.000000
#[4,]    0 0.0000000 0.2503962  1.0000000 -0.365899
#[5,]    0 0.0000000 0.0000000 -0.3658990  1.000000

# Third tau

res.fix5.3 <- fastGHS(x.sf.fix5,tau_sq = 0.125,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix5.3 <- cov2cor(res.fix5.3$theta)
theta.est.fix5.3[which(abs(theta.est.fix5.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.3!=0)
# 0.02828283
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.3!=0)
# 0.2214286
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.3!=0)
# 0.3131313
theta.est.fix5.3[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000  0.0000000  0.0000000
#[2,]    0 1.0000000 0.2640438  0.0000000  0.0000000
#[3,]    0 0.2640438 1.0000000  0.2569376  0.0000000
#[4,]    0 0.0000000 0.2569376  1.0000000 -0.3453257
#[5,]    0 0.0000000 0.0000000 -0.3453257  1.0000000

# As we see, if tau is of the same magntitude then the results are not too sensitive to the exact value - no significant problems arise. 

# EXAMPLE 6 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=500, p=100, larger partial correlations (0.183)

# Test sensitivity to choice of tau (within same magnitude), with stronger signal, small noise
n.fix6=500
p.fix6=100
set.seed(12345)
data.sf.fix6 = huge::huge.generator(n=n.fix6, d=p.fix6,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix6 = data.sf.fix6$theta # True adjacency matrix
theta.true.fix6 = data.sf.fix6$omega # The precision matrix
theta.true.fix6[which(theta.true.fix6<10e-5,arr.ind=T)]=0  
g.sf.fix6=graph.adjacency(data.sf.fix6$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix6 = data.sf.fix6$data # Observed attributes. nxp matrix.
x.sf.scaled.fix6= scale(x.sf.fix6) # Scale columns/variables.
s.sf.scaled.fix6 = cov(x.sf.scaled.fix6) # Empirical covariance matrix
data.sf.fix6$sparsity # True sparsity: 0.02
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix6[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1918105 0.0000000 0.0000000 0.0000000
#[2,] 0.1918105 1.0000000 0.1918105 0.0000000 0.0000000
#[3,] 0.0000000 0.1918105 1.0000000 0.1918105 0.1918105
#[4,] 0.0000000 0.0000000 0.1918105 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1918105 0.0000000 1.0000000

# First tau

res.fix6 <- fastGHS(x.sf.fix6,tau_sq = 0.7,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix6 <- cov2cor(res.fix6$theta)
theta.est.fix6[which(abs(theta.est.fix6) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix6!=0)
# 0.009090909
tailoredGlasso::precision(as.matrix(theta.true.fix6!=0), theta.est.fix6!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.fix6!=0), theta.est.fix6!=0)
# 0.4545455
theta.est.fix6[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2025924 0.0000000 0.0000000 0.0000000
#[2,] 0.2025924 1.0000000 0.2362132 0.0000000 0.0000000
#[3,] 0.0000000 0.2362132 1.0000000 0.2270435 0.1941258
#[4,] 0.0000000 0.0000000 0.2270435 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1941258 0.0000000 1.0000000

# Second tau

res.fix6.2 <- fastGHS(x.sf.fix6,tau_sq = 0.1,epsilon = 1e-5, fix_tau=TRUE)
theta.est.fix6.2 <- cov2cor(res.fix6.2$theta)
theta.est.fix6.2[which(abs(theta.est.fix6.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix6.2!=0)
# 0.007474747
tailoredGlasso::precision(as.matrix(theta.true.fix6!=0), theta.est.fix6.2!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.fix6!=0), theta.est.fix6.2!=0)
# 0.3737374
theta.est.fix6.2[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2021462 0.0000000 0.000000 0.0000000
#[2,] 0.2021462 1.0000000 0.2360934 0.000000 0.0000000
#[3,] 0.0000000 0.2360934 1.0000000 0.228576 0.1938736
#[4,] 0.0000000 0.0000000 0.2285760 1.000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1938736 0.000000 1.0000000


# With this strong signal, the choice of tau is much less important (but larger is better, as more edges are allowed in the estimated graph)



# EXAMPLE 7 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=50, p=100, larger partial correlations (0.192)

# Test sensitivity to initialisation

n.fix7=100
p.fix7=150
set.seed(12345)
data.sf.fix7 = huge::huge.generator(n=n.fix7, d=p.fix7,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix7 = data.sf.fix7$theta # True adjacency matrix
theta.true.fix7 = data.sf.fix7$omega # The precision matrix
theta.true.fix7[which(theta.true.fix7<10e-5,arr.ind=T)]=0  
g.sf.fix7=huge::graph.adjacency(data.sf.fix7$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix7 = data.sf.fix7$data # Observed attributes. nxp matrix.
x.sf.scaled.fix7= scale(x.sf.fix7) # Scale columns/variables.
s.sf.scaled.fix7 = cov(x.sf.scaled.fix7) # Empirical covariance matrix
data.sf.fix7$sparsity # True sparsity: 0.01333
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix7[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1826553 0.0000000 0.0000000 0.0000000
#[2,] 0.1826553 1.0000000 0.1826553 0.0000000 0.0000000
#[3,] 0.0000000 0.1826553 1.0000000 0.1826553 0.1826553
#[4,] 0.0000000 0.0000000 0.1826553 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1826553 0.0000000 1.0000000

# First initialisation (none)

res.fix7 <- fastGHS(x.sf.fix7,tau_sq = 0.05,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix7 <- cov2cor(res.fix7$theta)
theta.est.fix7[which(abs(theta.est.fix7) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix7!=0)
# 0.01243848
tailoredGlasso::precision(as.matrix(theta.true.fix7!=0), theta.est.fix7!=0)
# 0.3597122
tailoredGlasso::recall(as.matrix(theta.true.fix7!=0), theta.est.fix7!=0)
# 0.3355705
theta.est.fix7[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.2388029 0.0000000 0.0000000  0.0000000
#[2,] 0.2388029  1.0000000 0.1785712 0.0000000 -0.2205877
#[3,] 0.0000000  0.1785712 1.0000000 0.2034902  0.0000000
#[4,] 0.0000000  0.0000000 0.2034902 1.0000000  0.0000000
#[5,] 0.0000000 -0.2205877 0.0000000 0.0000000  1.0000000

# Second initialisation (random)

theta_init = matrix(runif(p.fix7^2,0,0.2),nrow=p.fix7)
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res.fix7.2 <- fastGHS(x.sf.fix7,theta=theta_init,tau_sq = 0.05,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix7.2 <- cov2cor(res.fix7.2$theta)
theta.est.fix7.2[which(abs(theta.est.fix7.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix7.2!=0)
# 0.01279642
tailoredGlasso::precision(as.matrix(theta.true.fix7!=0), theta.est.fix7.2!=0)
# 0.3566434
tailoredGlasso::recall(as.matrix(theta.true.fix7!=0), theta.est.fix7.2!=0)
# 0.3422819
theta.est.fix7.2[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.2388139 0.0000000 0.0000000  0.0000000
#[2,] 0.2388139  1.0000000 0.1786496 0.0000000 -0.2206028
#[3,] 0.0000000  0.1786496 1.0000000 0.2034793  0.0000000
#[4,] 0.0000000  0.0000000 0.2034793 1.0000000  0.0000000
#[5,] 0.0000000 -0.2206028 0.0000000 0.0000000  1.0000000

# Very similar results as to those with no initialisation!


# Third initialisation (random)

theta_init = matrix(runif(p.fix7^2,0,0.1),nrow=p.fix7)
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res.fix7.3 <- fastGHS(x.sf.fix7,theta=theta_init,tau_sq = 0.05,epsilon = 1e-3, fix_tau=TRUE)
theta.est.fix7.3 <- cov2cor(res.fix7.3$theta)
theta.est.fix7.3[which(abs(theta.est.fix7.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix7.3!=0)
# 0.01261745
tailoredGlasso::precision(as.matrix(theta.true.fix7!=0), theta.est.fix7.3!=0)
# 0.3687943
tailoredGlasso::recall(as.matrix(theta.true.fix7!=0), theta.est.fix7.3!=0)
# 0.3489933
theta.est.fix7.3[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000  0.0000000  0.0000000
#[2,]    0 1.0000000 0.2640438  0.0000000  0.0000000
#[3,]    0 0.2640438 1.0000000  0.2569376  0.0000000
#[4,]    0 0.0000000 0.2569376  1.0000000 -0.3453257
#[5,]    0 0.0000000 0.0000000 -0.3453257  1.0000000

# Very similar resultsto the others, only a bit better. 
# Sparsity is preserved in all cases.


=======
# With this stronger signal, the choice of tau is much less important (but larger is better, as more edges are allowed in the estimated graph)

