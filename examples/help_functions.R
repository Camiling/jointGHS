mutate.graph.2 = function(graph,fraction){
  # Swap a certain fraction of nodes in order to change graph
  # graph:      the huge.generate() object to mutate
  # fraction:   the fraction of edges to change
  prec.mat = graph$omega
  # Avoid rounding errors, which can happen with huge.generator
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  cov.mat = graph$sigma
  adj.mat = as.matrix(graph$theta)+0
  data=graph$data
  p = ncol(graph$omega)
  # If no nodes are to be swapped, simply return graph as is
  if(fraction==0){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    ans$data = data
    return(ans)
  }
  # We basically 'swap pairs of nodes' by switching their cols and rows
  # Indices of edge pairs
  edges = which(adj.mat==1,arr.ind=T) 
  # Avoid doubling it up
  edges=edges[1:(nrow(edges)/2),] 
  n.mutations = floor(nrow(edges)*fraction)
  # Sample nodes to give new edges
  nodes.add = sample(1:p,n.mutations) 
  # Sample nodes to remove edges from
  nodes.remove = edges[sample(1:nrow(edges),n.mutations),1] 
  # Must swap edges fpr adjacency matrix, precision matrix, covariance matrix and data matrix. 
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    # swap precision matrix rows. Then cols, the order does not matter
    prec.mat[nodes.remove[i],] = tmp.prec[nodes.add[i],]
    prec.mat[nodes.add[i],] = tmp.prec[nodes.remove[i],]
    tmp.prec=prec.mat
    # Swap precision matrix cols
    prec.mat[,nodes.remove[i]] = tmp.prec[,nodes.add[i]]
    prec.mat[,nodes.add[i]] = tmp.prec[,nodes.remove[i]]
    # Swap adjacency matrix rows
    adj.mat[nodes.remove[i],] = tmp.adj[nodes.add[i],]
    adj.mat[nodes.add[i],] = tmp.adj[nodes.remove[i],]
    tmp.adj = adj.mat
    # Swap adjacency matrix cols
    adj.mat[,nodes.remove[i]] = tmp.adj[,nodes.add[i]]
    adj.mat[,nodes.add[i]] = tmp.adj[,nodes.remove[i]]
    # Swap covariance mat rows
    cov.mat[nodes.remove[i],] = tmp.cov.mat[nodes.add[i],]
    cov.mat[nodes.add[i],] = tmp.cov.mat[nodes.remove[i],]
    tmp.cov.mat = cov.mat
    # Swap covariance mat cols
    cov.mat[,nodes.remove[i]] = tmp.cov.mat[,nodes.add[i]]
    cov.mat[,nodes.add[i]] = tmp.cov.mat[,nodes.remove[i]]
    # Swap data matrix cols (only the cols since it is nxp)
    data[,nodes.remove[i]] = tmp.dat[,nodes.add[i]]
    data[,nodes.add[i]] = tmp.dat[,nodes.remove[i]]
  }
  ans = list()
  ans$prec.mat = prec.mat
  ans$adj.mat = adj.mat
  ans$data = data
  return(ans)
}

mutate.graph= function(graph,fraction){
  # Mutate a given fraction of the edges of a graph. 
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change. 
  # We basically 'swap pairs of nodes' by switching their cols and rows. 
  prec.mat = graph$omega
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  cov.mat = graph$sigma
  adj.mat = graph$theta
  data=graph$data
  p = ncol(graph$omega)
  
  adj.mat.upper = adj.mat
  adj.mat.upper[lower.tri(adj.mat.upper)]=0
  diag(adj.mat.upper) =0
  edges = which(adj.mat.upper==1,arr.ind=T) # Edge pairs.
  n.mutations = floor(nrow(edges)*fraction)
  
  if(n.mutations==0 | is.na(n.mutations)){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    ans$data = data
    return(ans)
  }
  
  edges.to.change.ind = sample(1:nrow(edges),n.mutations) # We let the first index stay, then change the second one. 
  edges.to.change = edges[edges.to.change.ind,] # n.mutations x 2
  nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  nodes.remove = edges[sample(1:nrow(edges),n.mutations),1] # The nodes to 'leave out'
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    id.stay = edges.to.change[i,1]
    id.remove = edges.to.change[i,2]
    id.add=nodes.add[i]
    # Swap prec mat elements in rows. Then cols, the order does not matter!
    prec.mat[id.stay,id.add] = tmp.prec[id.stay,id.remove]
    prec.mat[id.stay,id.remove] = tmp.prec[id.stay,id.add]
    prec.mat[id.add,id.stay] = tmp.prec[id.remove,id.stay]
    prec.mat[id.remove,id.stay] = tmp.prec[id.add,id.stay]
    # swap adj mat rows
    adj.mat[id.stay,id.add] = tmp.adj[id.stay,id.remove]
    adj.mat[id.stay,id.remove] = tmp.adj[id.stay,id.add]
    adj.mat[id.add,id.stay] = tmp.adj[id.remove,id.stay]
    adj.mat[id.remove,id.stay] = tmp.adj[id.add,id.stay]
  }
  ans = list()
  ans$cov.mat=solve(prec.mat)
  ans$cov.mat[which(abs(ans$cov.mat)<1e-4,arr.ind = T)] = 0
  ans$prec.mat = prec.mat
  ans$adj.mat = adj.mat
  # Generate new data
  ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
  return(ans)
}

