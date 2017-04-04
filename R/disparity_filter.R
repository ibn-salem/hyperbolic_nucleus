
#Given an undirected network, it adds two properties to the edges. One containing the p-value for one end and the other for the other
disparity_filter <- function(net, deg.one.pval = 1){
  
  N <- vcount(net)
  
  E(net)$head.pval <- -1
  E(net)$tail.pval <- -1
  
  for(i in 1:N){
    #Get edges incident to node i
    edg <- incident(net, v = i)
    
    #Compute the node's degree
    k <- length(edg)
    
    if(k > 1){
      #Normalise the edge weights with respect to node i
      w <- edg$weight/sum(edg$weight)
      
      #Compute p-values for each edge based on beta distribution with shape parameters 1 and k-1
      pvals <- 1 - pbeta(w, shape1 = 1, shape2 = length(edg) - 1)
      
      #Identify already set pvals and put the computed one in the appropriate list
      #This depends on whether the current node is considered a head or a tail
      idx <- E(net)[edg]$head.pval < 0
      E(net)[edg[idx]]$head.pval <- pvals[idx]
      E(net)[edg[!idx]]$tail.pval <- pvals[!idx]
    }else{
      #Nodes of degree 1 get a p-value of deg.one.pval
      #Identify already set pvals and put the computed one in the appropriate list
      idx <- E(net)[edg]$head.pval < 0
      E(net)[edg[idx]]$head.pval <- deg.one.pval
      E(net)[edg[!idx]]$tail.pval <- deg.one.pval
    }
  }
  
  net$disparity.pval <- apply(cbind(E(net)$head.pval, E(net)$tail.pval), 1, min)
  
  return(net)
}

#After the application of the disparity filter, applies different pvalue thresholds and computes the fraction of remaining nodes and edges
#It does the same for a filter based on weights
analyse_filter <- function(bb, step = 0.01){
  disparity.cuts <- seq(max(E(bb)$disparity.pval), min(E(bb)$disparity.pval), -step)
  global.cuts <- seq(min(E(bb)$weight), max(E(bb)$weight), length.out = length(disparity.cuts))
  remaining <- data.frame(threshold = c(disparity.cuts, global.cuts), 
                          N = numeric(length = 2*length(disparity.cuts)), L = numeric(length = 2*length(disparity.cuts)), 
                          W = numeric(length = 2*length(disparity.cuts)), LCC = numeric(length = 2*length(disparity.cuts)),
                          filter = factor(rep(c("Disparity", "Global"), each = length(disparity.cuts)), levels = c("Disparity", "Global"), ordered = T))
  
  for(i in 1:nrow(remaining)){
    if(remaining$filter[i] == "Disparity"){
      g <- delete_edges(bb, edges = E(bb)[E(bb)$disparity.pval > remaining$threshold[i]])
    }else{
      g <- delete_edges(bb, edges = E(bb)[E(bb)$weight < remaining$threshold[i]])
    }
    g <- delete_vertices(g, v = which(degree(g) < 1))
    
    remaining$N[i] <- vcount(g)/vcount(bb)
    remaining$L[i] <- ecount(g)/ecount(bb)
    remaining$W[i] <- sum(E(g)$weight)/sum(E(bb)$weight)
    remaining$LCC[i] <- max(clusters(g)$csize)/vcount(bb)
  }
  
  return(remaining)
}

#Euclidean distances between a point x1 and a list of points x2
euc_dist <- function(x1, x2){
  return(sqrt((x1[1] - x2[, 1]) ^ 2 + (x1[2] - x2[, 2])^2))
}

#Find the threshold that provides the best trade-off between remaining fraction of nodes and edges
recommend_threshold <- function(bb, return.network = F){
  #Apply the filter over different thresholds
  remaining <- analyse_filter(bb, step = 0.01)
  
  #Compute distance between every point in the Lbb/Ltot<->Nbb/Ntot curve and the point (1, 0)
  d <- euc_dist(c(1, 0), remaining[remaining$filter == "Disparity", c("L", "N")])
  
  #Identify threshold
  threshold.info <- remaining[remaining$filter == "Disparity", ][which.max(d), ]
  
  #If requested, return filtered network and threshold. If not, return identified threshold
  if(return.network){
    filtered.net <- delete_edges(bb, edges = E(bb)[E(bb)$disparity.pval > threshold.info$threshold])
    return(list(filtered.net = filtered.net, threshold = threshold.info))
  }else{
    return(threshold.info)
  }
}

