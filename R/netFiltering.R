require(igraph) # for easy network generation and handling
require(cowplot) # to generate informative plots

source("R/disparity_filter.R")

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

#Construct network with promoter-promoter contact data
load("results/edge_lists.RData")
load("results/nodeDF.RData")

inData$weight <- inData$raw_count

#Construct network with all data
net <- graph_from_data_frame(d = inData, directed = F, vertices = nodeDF)

#Simplify (get rid of loops and multiple edges) and then take the largest connected component only
net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "sum")

#Apply Disparity Analysis [Serrano, Boguna & Vespignani (2009) PNAS 106(16)]

#Measure each promoter's disparity and plot degree vs disparity
disp <- get_node_disparity(net = net)
p.node.disp <- plot_degree_vs_disparity(net = net, node.disp = disp)

#Compute edge pvalues according to the disparity filter
edge.pvals <- get_edge_disparity_pvals(net = net)

#Analysis of different filters applied to the network
analysis <- analyse_disparity_filter(net = net, disparity.pval = edge.pvals, breaks = 50)

#Plot the analysis results
p.disp.filter <- plot_disparity_filter_analysis(disp.analysis = analysis)

#Put all plots together
disparity.plots <- plot_grid(p.node.disp, p.disp.filter$LvsN, p.disp.filter$WvsN, p.disp.filter$AvsLCC.tot, nrow = 2, ncol = 2, labels = letters[1:4])
save_plot("results/disparity_filter.pdf", disparity.plots, nrow = 2, ncol = 2, base_aspect_ratio = 1.3)

# #Filter network according to the recommended disparity threshold and get rid of isolated nodes
# filtered.net <- delete_edges(net, edges = E(net)[edge.pvals > analysis$threshold[analysis$recommended]])
# filtered.net <- delete_vertices(filtered.net, v = which(degree(filtered.net) < 1))
# 
# p.deg <- plot_degree_distr(network = filtered.net, bins = 500)
# 
# p.all <- plot_grid(plotlist = p.disp.filter, nrow = 1, ncol = 3, labels = letters[1:3])
# save_plot("results/disp_analysis.pdf", p.all, nrow = 1, ncol = 3, base_aspect_ratio = 1.3)

# save(disp, p.node.disp, edge.pvals, analysis, p.disp.filter, disparity.plots, p.deg, net, filtered.net, file = "results/disp_analysis.RData")
save(disp, p.node.disp, edge.pvals, analysis, p.disp.filter, disparity.plots, net, file = "results/disp_analysis.RData")


