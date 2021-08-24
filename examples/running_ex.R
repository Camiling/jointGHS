## Testing fastGHS

# The first simple running example

# Note that in these examples, tau is not fixed. 

# ---------------------------------------------
# Example 1: simple example, compare to the graphical lasso

library(fastGHS)
library(tailoredGlasso) # utils, download from https://github.com/Camiling/tailoredGlasso
library(huge)
library(network)
library(GGally)
set.seed(2020)
g <- huge::huge.generator(n=50,d=100,graph = 'scale-free')
X <- scale(g$data)
res <- fastGHS(X, epsilon = 1e-7, savepath = T)
# select the 2% edges with the largest partial correlation to get the correct sparsity
theta.est <- cov2cor(res$theta)
theta.est.off.diag <- theta.est
diag(theta.est.off.diag) <- NA
theta.est[which(abs(theta.est) < quantile(abs(theta.est.off.diag), 0.98,na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag), 0.98,na.rm = T)
#          98% 
# 4.870987e-14 
tailoredGlasso::sparsity(theta.est!=0)
# 0.02
g$sparsity
# 0.02
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est!=0)
# 0.3333333
tailoredGlasso::recall(as.matrix(g$theta!=0), theta.est!=0)
# 0.3333333
# Compare to graphical lasso
huge.g = huge(X,method='glasso', lambda = 0.35737)
huge.g$sparsity
# 0.02
precision(as.matrix(g$theta!=0),huge.g$icov[[1]]!=0)
# 0.3333333
recall(as.matrix(g$theta!=0), huge.g$icov[[1]]!=0)
# 0.3333333

# Same precision and recall for both methods. 

# -----------------------------------------
# Example 2: check that objective function increases with each step

res2 <- fastGHS(X, save_Q = T, epsilon = 1e-7)
res2$Q_vals

# It does

# -----------------------------------------
# Example 3: try to let in run longer with a very small epsilon

res3 <- fastGHS(X, savepath=T, save_Q = T, epsilon = 1e-300, maxitr = 100)
theta.est3 <- cov2cor(res3$theta)
theta.est.off.diag3 <- theta.est3
diag(theta.est.off.diag3) = NA
theta.est3[which(abs(theta.est3) < quantile(abs(theta.est.off.diag3), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::precision(as.matrix(g$theta!=0), theta.est3!=0)
# Same precision as for fewer runs

# Notably, elements still shrink strongly towards zero and are of very small magnitude

res3$tau_sq # Very small 
res3$Lambda_sq # Large diagonal elements, small off-diag elements


# ---------------------------------------------------
# Example 4: with groups that are unrelated -> should be separated

set.seed(2021)
g <- huge.generator(n=200,d=50,graph = 'scale-free') # Graph of first group
g2 <- huge.generator(n=200,d=50,graph = 'scale-free') # Graph of second group
g.full <- rbind(cbind(g$theta, matrix(rep(0,50^2), nrow=50)), cbind(matrix(rep(0,50^2), nrow=50),g2$theta)) # Combined 'true' graph with two connected components
sparsity(g.full!=0) # Sparsity of combined graph
X <- scale(cbind(g$data, g2$data)) # Combined data set
res4 <- fastGHS(X, group=c(rep(0,50), rep(1,50)), savepath=T, epsilon = 1e-5)
# Get network with correct sparsity
theta.est4 <- cov2cor(res4$theta)
theta.est.off.diag4 <- theta.est4
diag(theta.est.off.diag4) <- NA
theta.est4[which(abs(theta.est4) < quantile(abs(theta.est.off.diag4), 1-sparsity(g.full!=0),na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est4!=0)
gg <- theta.est4!=0

# Plot resulting fastGHS graph with groups marked by color
set.seed(1234)
net <- network::network(gg)
GGally::ggnet2(net,alpha=0.9,mode = "fruchtermanreingold",color = c(rep('deepskyblue2',50),rep('limegreen',50)))


# -----------------------------------------
# Example 5: test ICM

set.seed(2020) # Same graph in first example
g5 <- huge.generator(n=50,d=100,graph = 'scale-free')
X5 <- scale(g5$data)
res5 <- fastGHS(X5, savepath=F, epsilon = 1e-7, method='ICM')
theta.est5 <- cov2cor(res5$theta)
theta.est.off.diag5 <- theta.est5
diag(theta.est.off.diag5) <- NA
theta.est5[which(abs(theta.est5) < quantile(abs(theta.est.off.diag5), 0.98,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est5!=0)
# 0.02
g5$sparsity
# 0.02
tailoredGlasso::precision(as.matrix(g5$theta!=0), theta.est5!=0)
# 0.3333333
tailoredGlasso::recall(as.matrix(g5$theta!=0), theta.est5!=0)
# 0.3333333

# No difference from ECM result

tailoredGlasso::confusion.matrix(theta.est!=0, theta.est5!=0) # Complete edge agreement

# -------------------------------------------------
# Example 6: larger network

set.seed(222)
g6 <- huge.generator(n=500,d=300,graph = 'scale-free')
X6 <- scale(g6$data)
res6 <- fastGHS(X6, epsilon = 1e-200)
# select the 2% edges with the largest partial correlation to get the correct sparsity
theta.est6 <- cov2cor(res6$theta)
theta.est.off.diag6 <- theta.est6
diag(theta.est.off.diag6) <- NA
theta.est6[which(abs(theta.est6) < quantile(abs(theta.est.off.diag6), 1-g6$sparsity,na.rm = T), arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est6!=0)
# 0.006666667
g6$sparsity
# 0.006666667
tailoredGlasso::precision(as.matrix(g6$theta!=0), theta.est6!=0)
# 0.6187291
tailoredGlasso::recall(as.matrix(g6$theta!=0), theta.est6!=0)
# 0.6187291
