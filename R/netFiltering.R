require(igraph) # for easy network generation and handling
require(cowplot) # to generate informative plots

source("disparity_filter.R")

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

#Construct network with promoter-promoter contact data
load("../results/edge_lists.RData")

#Construct network with all data
net <- graph_from_data_frame(d = inData, directed = F)

#Simplify (get rid of loops and multiple edges) and get rid of isolated nodes
net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "min")
net <- delete_vertices(net, v = which(degree(net) < 1))
E(net)$weight <- E(net)$raw_count

#Apply Disparity Analysis [Serrano, Boguna & Vespignani (2009) PNAS 106(16)]

#Measure each promoter's disparity and plot degree vs disparity
disp <- get_node_disparity(net = net)
p.node.disp <- plot_degree_vs_disparity(net = net, node.disp = disp)

#Compute edge pvalues according to the disparity filter
edge.pvals <- get_edge_disparity_pvals(net = net)

#Analysis of different filters applied to the network
analysis <- analyse_disparity_filter(net = net, disparity.pval = edge.pvals, step = 0.01)

#Plot the analysis results
p.disp.filter <- plot_disparity_filter_analysis(disp.analysis = analysis)

#Put all plots together
disparity.plots <- plot_grid(p.node.disp, p.disp.filter$LvsN, p.disp.filter$WvsN, p.disp.filter$AvsLCC, nrow = 2, ncol = 2, labels = letters[1:4])
save_plot("../results/disparity_filter.pdf", disparity.plots, nrow = 2, ncol = 2, base_aspect_ratio = 1.3)

#Filter network according to the recommended disparity threshold and get rid of isolated nodes
filtered.net <- delete_edges(net, edges = E(net)[edge.pvals > analysis$threshold[analysis$recommended]])
filtered.net <- delete_vertices(filtered.net, v = which(degree(filtered.net) < 1))

p.deg <- plot_degree_distr(network = filtered.net)

p.all <- plot_grid(plotlist = p.disp.filter, nrow = 1, ncol = 3, labels = letters[1:3])
save_plot("../results/disp_analysis.pdf", p.all, nrow = 1, ncol = 3, base_aspect_ratio = 1.3)

save(disp, p.node.disp, edge.pvals, analysis, p.disp.filter, disparity.plots, p.deg, net, filtered.net, file = "../results/disp_analysis.RData")


