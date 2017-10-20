
require(igraph) # for easy network generation and handling
require(cowplot) # to generate informative plots

#' Generate a log10-linearly spaced sequence
#' 
#' Generates a sequence of \code{n} logarithmically spaced points between \code{s} and \code{e}, inclusive.
#' 
#' @param s numeric; The sequence starting point.
#' @param e numeric; The sequence ending point.
#' @param n integer; The number of points to be generated.
#' 
#' @return A numeric vector with \code{n} logarithmically spaced points between \code{s} and \code{e}.
#' 
log_seq <- function(s, e, n){
  s <- log10(s)
  e <- log10(e)
  return(exp(log(10) * seq(s, e, length.out = n)))
}

#' Plots the degree distribution of a network in the log-log space.
#' 
#' @param network An igraph object storing the network whose degree distribution is to be plotted.
#' @return An ggplot object
#' 
plot_degree_distr <- function(network, bins = 100){
  # Calculate degrees for each node
  d <- degree(network, mode = "all")
  
  # Bin the degrees using logarithmic spacing
  breaks = log_seq(min(d), max(d), bins + 1)
  cuts <- cut(d, breaks = breaks, labels = breaks[-length(breaks)], include.lowest = T)
  stats <- data.frame(deg = as.numeric(levels(cuts)), 
                      prob = as.numeric(table(cuts))/sum(as.numeric(table(cuts))), 
                      stringsAsFactors = F)
  stats <- stats[stats$prob > 0, ]
  
  return(ggplot(stats, aes_(~deg, ~prob)) + geom_point() + 
           scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4), 
                         labels = scales::trans_format("log10", scales::math_format())) + 
           scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4), 
                         labels = scales::trans_format("log10", scales::math_format())) + 
           annotation_logticks() + labs(x = "Node degree", y = "Probability") + theme_bw() +
           theme(panel.grid.minor = element_blank()))
}

#' Compares the clustering of a given network with that of 'epochs' ER graphs with the same number of nodes and edges
#' 
#' @param network An igraph object storing the network of interest.
#' @param epochs The number of ER graphs to consider in the comparison.
#' @param label A label for the given network in the resulting plot.
#' @return An ggplot object comparing the clustering coefficient of the given network and the average of 'epochs' ER graphs
#' 
compare_clustering_to_er <- function(network, epochs = 100, label = "Given network"){
  er.clust <- numeric(length = epochs)
  N <- vcount(network)
  E <- vcount(network)
  for(i in 1:epochs){
    g <- sample_gnm(N, E, directed = F, loops = F)
    er.clust[i] <- transitivity(g, type = "average")
  }
  df <- data.frame(net = factor(c(label, "ER graphs"), levels = c(label, "ER graphs"), ordered = T), 
                   clust = c(transitivity(network, type = "average"), mean(er.clust)),
                   err = c(0, sd(er.clust)))
  dodge <- position_dodge(width = 0.9)
  return(ggplot(df, aes(net, clust)) + geom_bar(position = dodge, stat = "identity") + geom_errorbar(aes(ymin = clust-err, ymax = clust+err), position = dodge, width = 0.25) + labs(x = "", y = "Clustering coefficient") + theme_bw())
}

#-------------------------------------------------------------------
# Some preliminary network analysis
#-------------------------------------------------------------------

load("results/edge_lists.RData")
load("results/noteDF.RData")

plots <- list()

inData$weight <- inData$raw_count

#Construct network with all data
net <- graph_from_data_frame(d = inData, directed = F, vertices = noteDF)

#Simplify (get rid of loops and multiple edges) and then take the largest connected component only
net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "sum")
co <- components(net)
# net <- induced_subgraph(net, vids = which(co$membership == which.max(co$csize)))

plots[[1]] <- plot_degree_distr(network = net, bins = 500)
plots[[2]] <- compare_clustering_to_er(network = net, epochs = 100, label = "Full Hi-C")

plots.full <- plot_grid(plotlist = plots, nrow = 1, ncol = 2, labels = letters[1:2])
save_plot("results/full_netAnalysis.pdf", plots.full, nrow = 1, ncol = 2, base_aspect_ratio = 1.3)

# #Construct network with filtered
# net <- graph_from_data_frame(d = fltDF, directed = F)
# 
# #Simplify (get rid of loops and multiple edges) and then take the largest connected component only
# net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "min")
# co <- clusters(net)
# net <- induced_subgraph(net, vids = which(co$membership == which.max(co$csize)))
# 
# plots[[1]] <- plot_degree_distr(network = net)
# plots[[2]] <- compare_clustering_to_er(network = net, epochs = 100, label = "Filtered Hi-C")
# 
# plots.filtered <- plot_grid(plotlist = plots, nrow = 1, ncol = 2, labels = letters[1:2])
# save_plot("results/filtered_netAnalysis.pdf", plots.filtered, nrow = 1, ncol = 2, base_aspect_ratio = 1.3)

