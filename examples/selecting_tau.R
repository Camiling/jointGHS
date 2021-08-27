library(fastGHS)
library(huge)
library(glasso)
library(igraph)
library(tailoredGlasso)

# Investigate how we can select tau by looking at how the sparsity plateaus

# EXAMPLE 1 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=200, p=50, larger partial correlations (0.229)

n.select=100
p.select=50
set.seed(12345)
data.sf.select = huge::huge.generator(n=n.select, d=p.select,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.select = data.sf.select$theta # True adjacency matrix
theta.true.select = data.sf.select$omega # The precision matrix
theta.true.select[which(theta.true.select<10e-5,arr.ind=T)]=0  
g.sf.select=graph.adjacency(data.sf.select$theta,mode="undirected",diag=F) # true igraph object
x.sf.select = data.sf.select$data # Observed attributes. nxp matrix.
data.sf.select$sparsity # True sparsity: 0.04
# Look at precision matrix (partial correlations)
cov2cor(theta.true.select[1:5,1:5])

taus = seq(1e-6,1,length.out = 100) # Actually tau square
res.select.ghs.list <- lapply(taus,FUN = function(l) fastGHS(x.sf.select,tau_sq = l,epsilon = 1e-3, fix_tau = T))
dev.new()
par(mfrow=c(2,2))
plot(taus, unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
plot(taus, unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
plot(taus, unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

dev.off()

# EXAMPLE 2 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=300, p=100, larger partial correlations (0.19)

n.select.2=300
p.select.2=100
set.seed(12345)
data.sf.select.2 = huge::huge.generator(n=n.select.2, d=p.select.2,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.select.2 = data.sf.select.2$theta # True adjacency matrix
theta.true.select.2 = data.sf.select.2$omega # The precision matrix
theta.true.select.2[which(theta.true.select.2<10e-5,arr.ind=T)]=0  
g.sf.select.2=graph.adjacency(data.sf.select.2$theta,mode="undirected",diag=F) # true igraph object
x.sf.select.2 = data.sf.select.2$data # Observed attributes. nxp matrix.
data.sf.select.2$sparsity # True sparsity: 0.02
# Look at precision matrix (partial correlations)
cov2cor(theta.true.select.2[1:5,1:5])

taus.2 = seq(1e-6,1,length.out = 100)
res.select.2.ghs.list <- lapply(taus.2,FUN = function(l) fastGHS(x.sf.select.2,tau_sq = l,epsilon = 1e-3, fix_tau = T))
dev.new()
par(mfrow=c(2,2))
plot(taus.2, unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
#par(new=T)
# As we see, the sparsity plateaus
plot(taus.2, unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.2!=0),abs(s$theta)>1e-5))), ,ylab='precision',type='l', col='blue')
#par(new=T)
plot(taus.2, unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.2!=0),abs(s$theta)>1e-5))), ,ylab='recall',type='l', col='red')

dev.off()

# EXAMPLE 3 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=200, p=150, larger partial correlations (0.18)

n.select.3=200
p.select.3=150
set.seed(12345)
data.sf.select.3 = huge::huge.generator(n=n.select.3, d=p.select.3,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.select.3 = data.sf.select.3$theta # True adjacency matrix
theta.true.select.3 = data.sf.select.3$omega # The precision matrix
theta.true.select.3[which(theta.true.select.3<10e-5,arr.ind=T)]=0  
g.sf.select.3=graph.adjacency(data.sf.select.3$theta,mode="undirected",diag=F) # true igraph object
x.sf.select.3 = data.sf.select.3$data # Observed attributes. nxp matrix.
data.sf.select.3$sparsity # True sparsity: 0.0133
# Look at precision matrix (partial correlations)
cov2cor(theta.true.select.3[1:5,1:5])

taus.3 = seq(1e-6,1,length.out = 50)
res.select.3.ghs.list <- lapply(taus.3,FUN = function(l) fastGHS(x.sf.select.3,tau_sq = l,epsilon = 1e-3, fix_tau = T))
dev.new()
par(mfrow=c(2,2))
plot(taus.3, unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
#par(new=T)
# As we see, the sparsity plateaus
plot(taus.3, unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.3!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
#par(new=T)
plot(taus.3, unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.3!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

dev.off()

# EXAMPLE 4 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=50, p=50, larger partial correlations (0.22)

n.select.4=50
p.select.4=50
set.seed(12345)
data.sf.select.4 = huge::huge.generator(n=n.select.4, d=p.select.4,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.select.4 = data.sf.select.4$theta # True adjacency matrix
theta.true.select.4 = data.sf.select.4$omega # The precision matrix
theta.true.select.4[which(theta.true.select.4<10e-5,arr.ind=T)]=0  
g.sf.select.4=graph.adjacency(data.sf.select.4$theta,mode="undirected",diag=F) # true igraph object
x.sf.select.4 = data.sf.select.4$data # Observed attributes. nxp matrix.
data.sf.select.4$sparsity # True sparsity: 0.04
# Look at precision matrix (partial correlations)
cov2cor(theta.true.select.4[1:5,1:5])

taus.4 = seq(1e-6,1,length.out = 100)
res.select.4.ghs.list <- lapply(taus.4,FUN = function(l) fastGHS(x.sf.select.4,tau_sq = l,epsilon = 1e-3, fix_tau = T))
dev.new()
par(mfrow=c(2,2))
plot(taus.4, unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
#par(new=T)
# As we see, the sparsity plateaus
plot(taus.4, unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.4!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
#par(new=T)
plot(taus.4, unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.4!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

dev.off()

# EXAMPLE 5 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=200, larger partial correlations (0.18)

n.select.5=100
p.select.5=200
set.seed(12345)
data.sf.select.5 = huge::huge.generator(n=n.select.5, d=p.select.5,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.select.5 = data.sf.select.5$theta # True adjacency matrix
theta.true.select.5 = data.sf.select.5$omega # The precision matrix
theta.true.select.5[which(theta.true.select.5<10e-5,arr.ind=T)]=0  
g.sf.select.5=graph.adjacency(data.sf.select.5$theta,mode="undirected",diag=F) # true igraph object
x.sf.select.5 = data.sf.select.5$data # Observed attributes. nxp matrix.
data.sf.select.5$sparsity # True sparsity: 0.01
# Look at precision matrix (partial correlations)
cov2cor(theta.true.select.5[1:5,1:5])

taus.5 = seq(1e-6,1,length.out = 50)
res.select.5.ghs.list <- lapply(taus.5,FUN = function(l) fastGHS(x.sf.select.5,tau_sq = l,epsilon = 1e-3, fix_tau = T))
dev.new()
par(mfrow=c(2,2))
plot(taus.5, unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
#par(new=T)
# As we see, the sparsity plateaus
plot(taus.5, unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.5!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
#par(new=T)
plot(taus.5, unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.5!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

dev.off()

# Plot all in one plot, to assess simultaneously ---------------------------------
pdf('~/Documents/Cambridge/PhD/fastGHS_files/fastGHS/examples/selectingTau.pdf',10,20)
par(mfrow=c(5,3))
plot(taus, unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
plot(taus, unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
plot(taus, unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

plot(taus.2, unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
plot(taus.2, unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.2!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
plot(taus.2, unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.2!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

plot(taus.3, unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
plot(taus.3, unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.3!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
plot(taus.3, unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.3!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

plot(taus.4, unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
plot(taus.4, unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.4!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
plot(taus.4, unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.4!=0),abs(s$theta)>1e-5))),ylab='recall',type='l', col='red')

plot(taus.5, unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5))), type='l',ylab='sparsity')
plot(taus.5, unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.5!=0),abs(s$theta)>1e-5))), ylab='precision',type='l', col='blue')
plot(taus.5, unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.5!=0),abs(s$theta)>1e-5))), ylab='recall',type='l', col='red')

dev.off()

# Test whether our proposed rule of thumb (smallest tau where the gradient of the sparsity is zero) works in these cases -------------------

spars.1 = unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5)))
spars.2 = unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5)))
spars.3 = unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5)))
spars.4 = unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5)))
spars.5 = unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::sparsity(abs(s$theta)>1e-5)))


# In real implementation: just save the two previous ones

get_tau_ind <- function(tau_vals, sparsities, eps=1e-3){
  n.tau = length(tau_vals)
  deriv = abs(sparsities[1:(n.tau-1)] - sparsities[2:n.tau])/abs(tau_vals[1:n.tau-1]- tau_vals[2:n.tau]) # Estimate the gradient 
  converged = deriv < eps
  converged_two_next = converged[1:(n.tau-2)] + converged[2:(n.tau-1)] + converged[3:n.tau]
  tau.ind = which(converged_two_next==3)[1]
  if (is.na(tau.ind)){
    return(get_tau_ind(tau_vals, sparsities, eps*10)) 
  }
  return(list(tau=tau_vals[tau.ind],tau_ind = tau.ind))
}

tau.ind.1 = get_tau_ind(taus, spars.1)
tau.ind.2 = get_tau_ind(taus.2, spars.2)
tau.ind.3 = get_tau_ind(taus.3, spars.3)
tau.ind.4 = get_tau_ind(taus.4, spars.4)
tau.ind.5 = get_tau_ind(taus.5, spars.5)

# Create plot with these points marked

prec.1 = unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select!=0),abs(s$theta)>1e-5)))
recall.1 = unlist(lapply(res.select.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select!=0),abs(s$theta)>1e-5)))
prec.2 = unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.2!=0),abs(s$theta)>1e-5)))
recall.2 = unlist(lapply(res.select.2.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.2!=0),abs(s$theta)>1e-5)))
prec.3 = unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.3!=0),abs(s$theta)>1e-5)))
recall.3 = unlist(lapply(res.select.3.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.3!=0),abs(s$theta)>1e-5)))
prec.4 = unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.4!=0),abs(s$theta)>1e-5)))
recall.4 = unlist(lapply(res.select.4.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.4!=0),abs(s$theta)>1e-5)))
prec.5 = unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::precision(as.matrix(theta.true.select.5!=0),abs(s$theta)>1e-5)))
recall.5 = unlist(lapply(res.select.5.ghs.list, FUN = function(s) tailoredGlasso::recall(as.matrix(theta.true.select.5!=0),abs(s$theta)>1e-5)))


pdf('~/Documents/Cambridge/PhD/fastGHS_files/fastGHS/examples/selectingTau_result.pdf',10,20)
par(mfrow=c(5,3))
plot(taus, spars.1, type='l',ylab='sparsity')
points(tau.ind.1$tau, spars.1[tau.ind.1$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus, prec.1, ylab='precision',type='l', col='blue')
points(tau.ind.1$tau, prec.1[tau.ind.1$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus, recall.1,ylab='recall',type='l', col='red')
points(tau.ind.1$tau, recall.1[tau.ind.1$tau_ind], col='darkorange', pch=16, cex=1.5)


plot(taus.2, spars.2, type='l',ylab='sparsity')
points(tau.ind.2$tau, spars.2[tau.ind.2$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.2, prec.2, ylab='precision',type='l', col='blue')
points(tau.ind.2$tau, prec.2[tau.ind.2$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.2, recall.2,ylab='recall',type='l', col='red')
points(tau.ind.2$tau, recall.2[tau.ind.2$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.3, spars.3, type='l',ylab='sparsity')
points(tau.ind.3$tau, spars.3[tau.ind.3$tau_ind], col='darkorange', pch=16)
plot(taus.3, prec.3, ylab='precision',type='l', col='blue')
points(tau.ind.3$tau, prec.3[tau.ind.3$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.3, recall.3, ylab='recall',type='l', col='red')
points(tau.ind.3$tau, recall.3[tau.ind.3$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.4, spars.4, type='l',ylab='sparsity')
points(tau.ind.4$tau, spars.4[tau.ind.4$tau_ind], col='darkorange', pch=16)
plot(taus.4, prec.4, ylab='precision',type='l', col='blue')
points(tau.ind.4$tau, prec.4[tau.ind.4$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.4, recall.4,ylab='recall',type='l', col='red')
points(tau.ind.4$tau, recall.4[tau.ind.4$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.5, spars.5, type='l',ylab='sparsity')
points(tau.ind.5$tau, spars.5[tau.ind.5$tau_ind], col='darkorange', pch=16)
plot(taus.5, prec.5, ylab='precision',type='l', col='blue')
points(tau.ind.5$tau, prec.5[tau.ind.5$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.5, recall.5,ylab='recall',type='l', col='red')
points(tau.ind.5$tau, recall.5[tau.ind.5$tau_ind], col='darkorange', pch=16, cex=1.5)

dev.off()

# Recreate plots with y axis going from 0 to 1 --------------------------------

pdf('~/Documents/Cambridge/PhD/fastGHS_files/fastGHS/examples/selectingTau_result_fullaxis.pdf',10,20)
par(mfrow=c(5,3))

plot(taus, spars.1, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select$sparsity+0.001,spars.1+0.001 )))
abline(data.sf.select$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.1$tau, spars.1[tau.ind.1$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus, prec.1, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.1$tau, prec.1[tau.ind.1$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus, recall.1,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.1$tau, recall.1[tau.ind.1$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.2, spars.2, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.2$sparsity+0.001,spars.2+0.001 )))
abline(data.sf.select.2$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.2$tau, spars.2[tau.ind.2$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.2, prec.2, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.2$tau, prec.2[tau.ind.2$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.2, recall.2,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.2$tau, recall.2[tau.ind.2$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.3, spars.3, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.3$sparsity+0.001,spars.3+0.001 )))
abline(data.sf.select.3$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.3$tau, spars.3[tau.ind.3$tau_ind], col='darkorange', pch=16)
plot(taus.3, prec.3, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.3$tau, prec.3[tau.ind.3$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.3, recall.3,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.3$tau, recall.3[tau.ind.3$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.4, spars.4, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.4$sparsity+0.001,spars.4+0.001 )))
abline(data.sf.select.4$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.4$tau, spars.4[tau.ind.4$tau_ind], col='darkorange', pch=16)
plot(taus.4, prec.4, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.4$tau, prec.4[tau.ind.4$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.4, recall.4, ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.4$tau, recall.4[tau.ind.4$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.5, spars.5, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.5$sparsity+0.001,spars.5+0.001 )))
abline(data.sf.select.5$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.5$tau, spars.5[tau.ind.5$tau_ind], col='darkorange', pch=16)
plot(taus.5, prec.5, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.5$tau, prec.5[tau.ind.5$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.5, recall.5,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.5$tau, recall.5[tau.ind.5$tau_ind], col='darkorange', pch=16, cex=1.5)

dev.off()


# Test alternative rule-of-thumb: tau corresponding to widest plateau (gradient zero) ---------------------------------------

get_tau_ind_max <- function(tau_vals, sparsities, eps=1e-3){
  n.tau = length(tau_vals)
  deriv = abs(sparsities[1:(n.tau-1)] - sparsities[2:n.tau])/abs(tau_vals[1:n.tau-1]- tau_vals[2:n.tau]) # Estimate the gradient 
  converged = deriv < eps
  converged_next = rep(0,length(converged))
  for(i in 1:length(converged)){
    j=0
    while(i+j <= length(converged) & converged[i+j]==T){
      converged_next[i] = converged_next[i] + 1
      j=j+1
    }
  }
  tau.ind = which.max(converged_next)
  if (is.na(tau.ind)){
    return(get_tau_ind_max(tau_vals, sparsities, eps*10)) 
  }
  return(list(tau=tau_vals[tau.ind],tau_ind = tau.ind))
}

tau.ind.1.max = get_tau_ind_max(taus, spars.1)
tau.ind.2.max = get_tau_ind_max(taus.2, spars.2)
tau.ind.3.max = get_tau_ind_max(taus.3, spars.3)
tau.ind.4.max = get_tau_ind_max(taus.4, spars.4)
tau.ind.5.max = get_tau_ind_max(taus.5, spars.5)

pdf('~/Documents/Cambridge/PhD/fastGHS_files/fastGHS/examples/selectingTau_max_result_fullaxis.pdf',10,20)
par(mfrow=c(5,3))
plot(taus, spars.1, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select$sparsity+0.001,spars.1+0.001 )))
abline(data.sf.select$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.1.max$tau, spars.1[tau.ind.1.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus, prec.1, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.1.max$tau, prec.1[tau.ind.1.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus, recall.1,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.1.max$tau, recall.1[tau.ind.1.max$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.2, spars.2, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.2$sparsity+0.001,spars.2+0.001 )))
abline(data.sf.select.2$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.2.max$tau, spars.2[tau.ind.2.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.2, prec.2, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.2.max$tau, prec.2[tau.ind.2.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.2, recall.2,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.2.max$tau, recall.2[tau.ind.2.max$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.3, spars.3, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.3$sparsity+0.001,spars.3+0.001 )))
abline(data.sf.select.3$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.3.max$tau, spars.3[tau.ind.3.max$tau_ind], col='darkorange', pch=16)
plot(taus.3, prec.3, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.3.max$tau, prec.3[tau.ind.3.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.3, recall.3,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.3.max$tau, recall.3[tau.ind.3.max$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.4, spars.4, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.4$sparsity+0.001,spars.4+0.001 )))
abline(data.sf.select.4$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.4.max$tau, spars.4[tau.ind.4.max$tau_ind], col='darkorange', pch=16)
plot(taus.4, prec.4, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.4.max$tau, prec.4[tau.ind.4.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.4, recall.4, ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.4.max$tau, recall.4[tau.ind.4.max$tau_ind], col='darkorange', pch=16, cex=1.5)

plot(taus.5, spars.5, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.5$sparsity+0.001,spars.5+0.001 )))
abline(data.sf.select.5$sparsity, 0, col='limegreen', lty=2)
points(tau.ind.5.max$tau, spars.5[tau.ind.5.max$tau_ind], col='darkorange', pch=16)
plot(taus.5, prec.5, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.ind.5.max$tau, prec.5[tau.ind.5.max$tau_ind], col='darkorange', pch=16, cex=1.5)
plot(taus.5, recall.5,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.ind.5.max$tau, recall.5[tau.ind.5.max$tau_ind], col='darkorange', pch=16, cex=1.5)

dev.off()


# Test Williams et al. approach to selecting tau (based on Piironen & Vehtari) -------------------------

# Prior expectation of sparsity
alpha.prior = seq(0,0.1,length.out = 100)
tau.prior.1 = alpha.prior/(1-alpha.prior)/sqrt(n.select)
cbind(alpha.prior,tau.prior.1)

tau.prior.2 = alpha.prior/(1-alpha.prior)/sqrt(n.select.2)
cbind(alpha.prior,tau.prior.2)

tau.prior.3 = alpha.prior/(1-alpha.prior)/sqrt(n.select.3)
cbind(alpha.prior,tau.prior.3)

tau.prior.4 = alpha.prior/(1-alpha.prior)/sqrt(n.select.4)
cbind(alpha.prior,tau.prior.4)

tau.prior.5 = alpha.prior/(1-alpha.prior)/sqrt(n.select.5)
cbind(alpha.prior,tau.prior.5)

# They are all very small....

# Recreate plots with these values marked out

tab.1 = as.data.frame(cbind(taus, spars.1))
ind.tau.1 = sapply(tau.prior.1, FUN=function(t) which.min(abs(t^2-taus)))

tab.2 = as.data.frame(cbind(taus.2, spars.2))
ind.tau.2 = sapply(tau.prior.2, FUN=function(t) which.min(abs(t^2-taus)))

tab.3 = as.data.frame(cbind(taus.3, spars.3))
ind.tau.3 = sapply(tau.prior.3, FUN=function(t) which.min(abs(t^2-taus)))

tab.4 = as.data.frame(cbind(taus.4, spars.4))
ind.tau.4 = sapply(tau.prior.4, FUN=function(t) which.min(abs(t^2-taus)))

tab.5 = as.data.frame(cbind(taus.5, spars.5))
ind.tau.5 = sapply(tau.prior.5, FUN=function(t) which.min(abs(t^2-taus)))

pdf('~/Documents/Cambridge/PhD/fastGHS_files/fastGHS/examples/selectingTau_Williams.pdf',10,20)
par(mfrow=c(5,3))
plot(taus, spars.1, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select$sparsity+0.001,spars.1+0.001 )))
abline(data.sf.select$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.1, spars.1[ind.tau.1],col='darkorchid2', pch=16, cex=1.5)
plot(taus, prec.1, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.1, prec.1[ind.tau.1],col='darkorchid2', pch=16, cex=1.5)
plot(taus, recall.1,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.1, recall.1[ind.tau.1],col='darkorchid2', pch=16, cex=1.5)

plot(taus.2, spars.2, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.2$sparsity+0.001,spars.2+0.001 )))
abline(data.sf.select.2$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.2, spars.2[ind.tau.2],col='darkorchid2', pch=16, cex=1.5)
plot(taus.2, prec.2, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.2, prec.2[ind.tau.2],col='darkorchid2', pch=16, cex=1.5)
plot(taus.2, recall.2,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.2, recall.2[ind.tau.2],col='darkorchid2', pch=16, cex=1.5)


plot(taus.3, spars.3, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.3$sparsity+0.001,spars.3+0.001 )))
abline(data.sf.select.3$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.3, spars.3[ind.tau.3],col='darkorchid2', pch=16, cex=1.5)
plot(taus.3, prec.3, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.3, prec.3[ind.tau.3],col='darkorchid2', pch=16, cex=1.5)
plot(taus.3, recall.3,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.3, recall.3[ind.tau.3],col='darkorchid2', pch=16, cex=1.5)


plot(taus.4, spars.4, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.4$sparsity+0.001,spars.4+0.001 )))
abline(data.sf.select.4$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.4, spars.4[ind.tau.4],col='darkorchid2', pch=16, cex=1.5)
plot(taus.4, prec.4, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.4, prec.4[ind.tau.4],col='darkorchid2', pch=16, cex=1.5)
plot(taus.4, recall.4, ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.4, recall.4[ind.tau.4],col='darkorchid2', pch=16, cex=1.5)

plot(taus.5, spars.5, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.5$sparsity+0.001,spars.5+0.001 )))
abline(data.sf.select.5$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.5, spars.5[ind.tau.5],col='darkorchid2', pch=16, cex=1.5)
plot(taus.5, prec.5, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.5, prec.5[ind.tau.5],col='darkorchid2', pch=16, cex=1.5)
plot(taus.5, recall.5,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.5, recall.5[ind.tau.5],col='darkorchid2', pch=16, cex=1.5)

dev.off()

# Also create plots of alpha vs tau to illustrate how it changes

plot(alpha.prior, tau.prior.1)
plot(alpha.prior,tau.prior.2)
plot(alpha.prior, tau.prior.3)
plot(alpha.prior,tau.prior.4)
plot(alpha.prior, tau.prior.5)

# Increases linearly?


# Recreate plots with tau selected with only alpha = true sparsity marked

alpha.prior.1 = data.sf.select$sparsity
tau.prior.1 = alpha.prior.1/(1-alpha.prior.1)/sqrt(n.select)
cbind(alpha.prior.1,tau.prior.1)
res.alpha.1 = fastGHS(x.sf.select,tau_sq = tau.prior.1^2 ,epsilon = 1e-3, fix_tau = T)

alpha.prior.2 = data.sf.select.2$sparsity
tau.prior.2 = alpha.prior.2/(1-alpha.prior.2)/sqrt(n.select.2)
cbind(alpha.prior.2,tau.prior.2)
res.alpha.2 = fastGHS(x.sf.select.2,tau_sq = tau.prior.2^2 ,epsilon = 1e-3, fix_tau = T)

alpha.prior.3 = data.sf.select.3$sparsity
tau.prior.3 = alpha.prior.3/(1-alpha.prior.3)/sqrt(n.select.3)
cbind(alpha.prior.3,tau.prior.3)
res.alpha.3 = fastGHS(x.sf.select.3,tau_sq = tau.prior.3^2 ,epsilon = 1e-3, fix_tau = T)

alpha.prior.4 = data.sf.select.4$sparsity
tau.prior.4 = alpha.prior.4/(1-alpha.prior.4)/sqrt(n.select.4)
cbind(alpha.prior.4,tau.prior.4)
res.alpha.4 = fastGHS(x.sf.select.4,tau_sq = tau.prior.4^2 ,epsilon = 1e-3, fix_tau = T)

alpha.prior.5 = data.sf.select.5$sparsity
tau.prior.5 = alpha.prior.5/(1-alpha.prior.5)/sqrt(n.select.5)
cbind(alpha.prior.5,tau.prior.5)
res.alpha.5 = fastGHS(x.sf.select.5,tau_sq = tau.prior.5^2 ,epsilon = 1e-3, fix_tau = T)

# As the tau values are so small, we must use a higher threshold in order to get readable results. 
hist(abs(res.alpha.1$theta[!diag(p.select)]))
# not a clear distinction between zero- and non-zero values: not a good tau value. 

pdf('~/Documents/Cambridge/PhD/fastGHS_files/fastGHS/examples/selectingTau_Williams_truesize.pdf',10,20)
par(mfrow=c(5,3))
plot(taus, spars.1, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select$sparsity+0.001,spars.1+0.001 )))
abline(data.sf.select$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.1, sparsity(abs(res.alpha.1$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus, prec.1, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.1, precision(theta.true.select!=0,abs(res.alpha.1$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus, recall.1,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.1, recall(theta.true.select!=0,abs(res.alpha.1$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)

plot(taus.2, spars.2, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.2$sparsity+0.001,spars.2+0.001 )))
abline(data.sf.select.2$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.2, sparsity(abs(res.alpha.2$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.2, prec.2, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.2, precision(theta.true.select.2!=0,abs(res.alpha.2$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.2, recall.2,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.2, recall(theta.true.select.2!=0,abs(res.alpha.2$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)


plot(taus.3, spars.3, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.3$sparsity+0.001,spars.3+0.001 )))
abline(data.sf.select.3$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.3, sparsity(abs(res.alpha.3$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.3, prec.3, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.3, precision(theta.true.select.3!=0,abs(res.alpha.3$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.3, recall.3,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.3, recall(theta.true.select.3!=0,abs(res.alpha.3$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)


plot(taus.4, spars.4, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.4$sparsity+0.001,spars.4+0.001 )))
abline(data.sf.select.4$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.4, sparsity(abs(res.alpha.4$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.4, prec.4, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.4, precision(theta.true.select.4!=0,abs(res.alpha.4$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.4, recall.4, ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.4, recall(theta.true.select.4!=0,abs(res.alpha.4$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)

plot(taus.5, spars.5, type='l',ylab='sparsity', ylim=c(0,max(data.sf.select.5$sparsity+0.001,spars.5+0.001 )))
abline(data.sf.select.5$sparsity, 0, col='limegreen', lty=2)
points(tau.prior.5, sparsity(abs(res.alpha.5$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.5, prec.5, ylab='precision',type='l', col='deepskyblue2', ylim=c(0,1))
points(tau.prior.5, precision(theta.true.select.5!=0,abs(res.alpha.5$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)
plot(taus.5, recall.5,ylab='recall',type='l', col='deeppink3', ylim=c(0,1))
points(tau.prior.5, recall(theta.true.select.5!=0,abs(res.alpha.5$theta)>1e-4),col='darkorchid2', pch=16, cex=1.5)

dev.off()





