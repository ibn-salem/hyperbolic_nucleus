
require(igraph) # for easy network generation and handling
require(cowplot) # to generate informative plots


#' Plots the degree distribution of a network in the log-log space.
#' 
#' @param network An igraph object storing the network whose degree distribution is to be plotted.
#' @return An ggplot object
#' 
plot_degree_distr <- function(network){
  # calculate degree
  d <- degree(network, mode = "all")
  dd <- degree_distribution(network, mode = "all", cumulative = FALSE)
  degs <- 1:max(d)
  prob <- dd[-1]
  # delete blank values
  nonzero.position = which(prob != 0)
  prob <- prob[nonzero.position]
  degs <- degs[nonzero.position]
  return(ggplot(data.frame(degs, prob), aes(degs, prob)) + geom_point() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks() + labs(x = "Node degree", y = "Probability") + theme_bw())
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

plots <- list()

#Construct network with all data
net <- graph_from_data_frame(d = inData, directed = F)

#Simplify (get rid of loops and multiple edges) and then take the largest connected component only
net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "min")
co <- clusters(net)
net <- induced_subgraph(net, vids = which(co$membership == which.max(co$csize)))

plots[[1]] <- plot_degree_distr(network = net)
plots[[2]] <- compare_clustering_to_er(network = net, epochs = 100, label = "Full Hi-C")

plots.full <- plot_grid(plotlist = plots, nrow = 1, ncol = 2, labels = letters[1:2])
save_plot("results/full_netAnalysis.pdf", plots.full, nrow = 1, ncol = 2, base_aspect_ratio = 1.3)

#Construct network with filtered
net <- graph_from_data_frame(d = fltDF, directed = F)

#Simplify (get rid of loops and multiple edges) and then take the largest connected component only
net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "min")
co <- clusters(net)
net <- induced_subgraph(net, vids = which(co$membership == which.max(co$csize)))

plots[[1]] <- plot_degree_distr(network = net)
plots[[2]] <- compare_clustering_to_er(network = net, epochs = 100, label = "Filtered Hi-C")

plots.filtered <- plot_grid(plotlist = plots, nrow = 1, ncol = 2, labels = letters[1:2])
save_plot("results/filtered_netAnalysis.pdf", plots.filtered, nrow = 1, ncol = 2, base_aspect_ratio = 1.3)

