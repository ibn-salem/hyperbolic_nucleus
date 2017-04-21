require(igraph)

#' Disparity measure of network nodes
#' 
#' Computes the disparity measure of each node in a given weighted network.
#' 
#' @param net igraph; A weighted undirected network.
#' 
#' @return A numeric vector of length \code{vcount(net)} with the disparity measure of each node.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone of complex weighted networks. \emph{PNAS} 106(16).
#' @references Barthelemy, M. et al. (2003) Spatial structure of the internet traffic. \emph{Physica A} 319.
#' 
#' @examples
#' # Compute the disparity of the nodes in the included US Airports network
#' disp <- get_node_disparity(net = air)
#' 
#' @export
#' @import igraph
#' 
get_node_disparity <- function(net){
  N <- vcount(net)
  disp <- numeric(length = N)
  for(i in 1:N){
    # Get edges incident to node i
    edg <- incident(net, v = i)
    
    # Compute the node's degree
    k <- length(edg)
    
    # Normalise the edge weights with respect to node i
    w <- edg$weight/sum(edg$weight)
    
    # Compute the disparity for node i
    disp[i] <- k*sum(w^2)
  }
  return(disp)
}

#' Disparity measures of the null model
#' 
#' Given a degree sequence, computes the average and standard deviation of the disparities for the null corresponding to each degree.
#' 
#' @param degrees numeric; A numeric vector containing a node degree sequence.
#' 
#' @return List with the following elements:
#' \item{mu}{Numeric vector with the null average disparities for the given degrees.}
#' \item{sd}{Numeric vector with the null standard deviation of disparities for the given degrees.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone of complex weighted networks. \emph{PNAS} 106(16).
#'
#' @examples
#' # Compute the null disparities for the spectrum of node degrees of the included US Airports network
#' null.disp <- get_node_null_disparity(degrees = seq(from = 1, to = max(degree(air)), by = 1))
#' 
#' @export
#' 
get_node_null_disparity <- function(degrees){
  mu <- 2*degrees/(degrees + 1)
  sig <- sqrt((degrees^2)*((20 + 4*degrees)/((degrees + 1)*(degrees + 2)*(degrees + 3)) - (4/((degrees + 1)^2))))
  return(list(mu = mu, sd = sig))
}

#' Visualise node disparities as a function of node degree
#' 
#' Given a  weighted undirected network and node disparities precomputed with function \code{get_node_disparity}, 
#' generates a ggplot comparing node disparities as a function of node degree.
#' Reference lines for null disparities, perfect node homogeneity and perfect node heterogeneity are also shown.
#' 
#' @param net igraph; A weighted undirected network.
#' @param node.disp numeric; Numeric vector of length \code{vcount(net)} with the disparity measure of each node in the network.
#' 
#' @return A ggplot comparing node disparities as a function of node degree.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone of complex weighted networks. \emph{PNAS} 106(16).
#'
#' @examples
#' # Compute the disparity of the nodes in the included US Airports network
#' disp <- get_node_disparity(net = air)
#' # Plot degree vs disparities
#' p <- plot_degree_vs_disparity(net = air, node.disp = disp)
#' 
#' @export
#' @import igraph
#' 
plot_degree_vs_disparity <- function(net, node.disp){
  null.seq <- seq(from = 1, to = max(degree(net)), by = 1)
  null.disp <- get_node_null_disparity(degrees = null.seq)
  
  df.main <- data.frame(deg = degree(net), disp = node.disp, stringsAsFactors = F)
  df.null <- data.frame(deg = null.seq, disp = null.disp$mu + 2*null.disp$sd, stringsAsFactors = F)
  
  p <- ggplot(df.main, aes(x = deg, y = disp)) + geom_point(colour = "#7fbf7b") +
    geom_hline(yintercept = 1, linetype = 2, colour = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "blue") +
    geom_line(data = df.null, aes(x = deg, y = disp), colour = "red") + 
    coord_cartesian(xlim = c(1, max(degree(net))), ylim = c(1, max(degree(net)))) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + annotation_logticks() + labs(x = "Node degree", y = "Node disparity") + theme_bw() +
    theme(panel.grid.minor = element_blank())
}

#' Disparity filter for a weighted undirected network
#' 
#' For each weighted edge in the given network, it computes two p-values representing its local significance for the two nodes it touches.
#' This function returns only the minium p-value for the two that are computed.
#' 
#' @param net igraph; A weighted undirected network.
#' @param deg.one.pval numeric; By default, edges of degree-one nodes are assigned a p-value of 1. This can be changed through this parameter.
#' 
#' @return A numeric vector of length \code{ecount(net)} with a p-value for each weighted edge of the input network.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone of complex weighted networks. \emph{PNAS} 106(16).
#' 
#' @examples
#' # Get disparity p-values for the edges of the included US Airports network
#' edge.pvals <- get_edge_disparity_pvals(net = air)
#' 
#' @export
#' @import igraph
#' 
get_edge_disparity_pvals <- function(net, deg.one.pval = 1){
  
  N <- vcount(net)
  
  E(net)$head.pval <- -1
  E(net)$tail.pval <- -1
  
  for(i in 1:N){
    # Get edges incident to node i
    edg <- incident(net, v = i)
    
    # Compute the node's degree
    k <- length(edg)
    
    if(k > 1){
      # Normalise the edge weights with respect to node i
      w <- edg$weight/sum(edg$weight)
      
      # Compute p-values for each edge based on beta distribution with shape parameters 1 and k-1
      pvals <- 1 - pbeta(w, shape1 = 1, shape2 = length(edg) - 1)
      
      # Identify already set pvals and put the computed one in the appropriate list
      # This depends on whether the current node is considered a head or a tail
      idx <- E(net)[edg]$head.pval < 0
      E(net)[edg[idx]]$head.pval <- pvals[idx]
      E(net)[edg[!idx]]$tail.pval <- pvals[!idx]
    }else{
      # Nodes of degree 1 get a p-value of deg.one.pval
      # Identify already set pvals and put the computed one in the appropriate list
      idx <- E(net)[edg]$head.pval < 0
      E(net)[edg[idx]]$head.pval <- deg.one.pval
      E(net)[edg[!idx]]$tail.pval <- deg.one.pval
    }
  }
  
  # Obtain the minimum p-value from the two computed for each edge
  disparity.pval <- apply(cbind(E(net)$head.pval, E(net)$tail.pval), 1, min)
  
  return(disparity.pval)
}

#' Topological analysis of the disparity filter
#' 
#' After disparity p-values have been obtained, applies different p-value thresholds to the given network and computes the fraction of 
#' remaining nodes, edges and size of the largest connected component from the perspective of the original network.
#' It does the same for a filter based on weights (global filter).
#' In addition, it recommends a disparity p-value to filter the network, keeping as many nodes as possible and the most significant edges.
#' 
#' @param net igraph; The weighted undirected network to which the disparity filter was applied.
#' @param disparity.pval numeric; The vector of edge p-values obtained with function \code{get_edge_disparity_pvals} for network \code{net}.
#' @param step numeric; The size of the step for the sequence of values between the maximum and minimum edge p-values.
#' 
#' @return Data frame with the following columns:
#' \item{threshold}{The different p-value thresholds considered in the analysis.}
#' \item{N}{The remaining fraction of nodes at each threshold.}
#' \item{L}{The remaining fraction of links at each threshold.}
#' \item{W}{The remaining fraction of the total weight at each threshold.}
#' \item{LCC}{Fraction of nodes from the original network in the largest connected component at each threshold.}
#' \item{recommended}{TRUE if the threshold is the recommended one to filter the network, FALSE otherwise.}
#' \item{filter}{Whether the entry of the data frame corresponds to the Disparity or the Global filter.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone of complex weighted networks. \emph{PNAS} 106(16).
#' @references García-Pérez, G. et al. (2016) The hidden hyperbolic geometry of international trade: World Trade Atlas 1870–2013. 
#' \emph{Scientific Reports} 6(33441).
#'
#' @examples
#' # Get disparity p-values for the edges of the included US Airports network
#' edge.pvals <- get_edge_disparity_pvals(net = air)
#' # Analyse the topology of the networks resulting from the application of different disparity filters
#' analysis <- analyse_disparity_filter(net = air, disparity.pval = edge.pvals, step = 0.01)
#' 
#' @export
#' @import igraph
#' 
analyse_disparity_filter <- function(net, disparity.pval, step = 0.01){
  disparity.cuts <- seq(from = 1, to = 0, by = -step)
  global.cuts <- seq(from = min(E(net)$weight), to = max(E(net)$weight), length.out = length(disparity.cuts))
  
  remaining <- data.frame(threshold = c(disparity.cuts, global.cuts), 
                          N = numeric(length = 2*length(disparity.cuts)), L = numeric(length = 2*length(disparity.cuts)), 
                          W = numeric(length = 2*length(disparity.cuts)), LCC = numeric(length = 2*length(disparity.cuts)),
                          recommended = logical(length = 2*length(disparity.cuts)), filter = factor(rep(c("Disparity", "Global"), 
                                                                                                        each = length(disparity.cuts)), 
                                                                                                    levels = c("Disparity", "Global"), ordered = T))
  
  for(i in 1:nrow(remaining)){
    if(remaining$filter[i] == "Disparity"){
      g <- delete_edges(net, edges = E(net)[disparity.pval > remaining$threshold[i]])
    }else{
      g <- delete_edges(net, edges = E(net)[E(net)$weight < remaining$threshold[i]])
    }
    g <- delete_vertices(g, v = which(degree(g) < 1))
    
    remaining$N[i] <- vcount(g)/vcount(net)
    remaining$L[i] <- ecount(g)/ecount(net)
    remaining$W[i] <- sum(E(g)$weight)/sum(E(net)$weight)
    remaining$LCC[i] <- ifelse(length(clusters(g)$csize) > 0, max(clusters(g)$csize)/vcount(net), 0)
  }
  
  # Compute the vertical distance between every point in the Lbb/Ltot<->Nbb/Ntot curve and the diagonal (where Lbb/Ltot == Nbb/Ntot)
  diag.vals <- remaining$L[remaining$filter == "Disparity"]
  d <- sqrt((remaining$L[remaining$filter == "Disparity"] - diag.vals)^2 + (remaining$N[remaining$filter == "Disparity"] - diag.vals)^2)
  # The Lbb/Ltot<->Nbb/Ntot pair that maximises the vertical distance to the diagonal is set as the recommended threshold
  remaining$recommended[which.max(d)] <- TRUE
  
  return(remaining)
}

#' Visual analysis of the disparity filter applied to a weighted undirected network
#' 
#' After the application of function \code{analyse_disparity_filter} to a weighted undirected network,
#' it generates three ggplots showing how the remaining fraction of nodes changes as a function of the
#' remaining fraction of links and total weight as more stringent filters are applied to the network.
#' The third plot shows how the fraction of nodes in the largest connected component of the network 
#' changes as a function of the different filters applied.
#' All three plots indicate the final situation of the network for the recommended disparity filter.
#' 
#' @param analysis data frame; The data frame resulting from the application of function \code{analyse_disparity_filter}.
#' 
#' @return A list of three ggplots:
#' \item{LvsN}{The remaining fraction of nodes as a function of the remaining fraction of links.}
#' \item{WvsN}{The remaining fraction of nodes as a function of the remaining fraction of total weight.}
#' \item{AvsLCC}{Fraction of nodes from the original network in the LCC as a function of different p-value thresholds.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone of complex weighted networks. \emph{PNAS} 106(16).
#' @references García-Pérez, G. et al. (2016) The hidden hyperbolic geometry of international trade: World Trade Atlas 1870–2013. 
#' \emph{Scientific Reports} 6(33441).
#'
#' @examples
#' # Get disparity p-values for the edges of the included US Airports network
#' edge.pvals <- get_edge_disparity_pvals(net = air)
#' # Analyse the topology of the networks resulting from the application of different disparity filters
#' analysis <- analyse_disparity_filter(net = air, disp.pvals = edge.pvals, step = 0.01)
#' # Plot the results of the analysis
#' p <- plot_disparity_filter_analysis(disp.analysis = analysis)
#' 
#' @export
#' @import ggplot2
#' 
plot_disparity_filter_analysis <- function(disp.analysis){
  p <- vector(mode = "list", length = 3)
  names(p) <- c("LvsN", "WvsN", "AvsLCC")
  
  p$LvsN <- ggplot(disp.analysis, aes(L, N, colour = filter, shape = filter)) + geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_abline(slope = -1, intercept = 0, linetype = 2, colour = "blue") +
    geom_vline(xintercept = disp.analysis$L[disp.analysis$recommended], linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp.analysis$N[disp.analysis$recommended], linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = expression(L[bb]/L[tot]), y = expression(N[bb]/N[tot])) + theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), legend.justification = c(0,0), legend.position = c(0,0))
  
  p$WvsN <- ggplot(disp.analysis, aes(W, N, colour = filter, shape = filter)) + geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_vline(xintercept = disp.analysis$W[disp.analysis$recommended], linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp.analysis$N[disp.analysis$recommended], linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = expression(W[bb]/W[tot]), y = expression(N[bb]/N[tot])) + theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), legend.justification = c(0,0), legend.position = c(0,0))
  
  #Normalise the weight for the Global filter
  disp.analysis$threshold[disp.analysis$filter == "Global"] <- (disp.analysis$threshold[disp.analysis$filter == "Global"] - 
                                                                  min(disp.analysis$threshold[disp.analysis$filter == "Global"])) /
    (max(disp.analysis$threshold[disp.analysis$filter == "Global"]) - 
       min(disp.analysis$threshold[disp.analysis$filter == "Global"]))
  
  p$AvsLCC <- ggplot(disp.analysis, aes(threshold, LCC, colour = filter, shape = filter)) + geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_vline(xintercept = disp.analysis$threshold[disp.analysis$recommended], linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp.analysis$LCC[disp.analysis$recommended], linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = "Disparity p-values/Normalised weight", y = expression(paste("Fraction of ", N[tot], " in ", LCC[bb]))) + theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), legend.justification = c(0,0), legend.position = c(0,0.5))
  
  return(p)
}
