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

plots <- list()

#Construct network with all data
net <- graph_from_data_frame(d = inData, directed = F)

#Simplify (get rid of loops and multiple edges) and then take the largest connected component only
net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb = "min")
co <- clusters(net)
net <- induced_subgraph(net, vids = which(co$membership == which.max(co$csize)))
E(net)$weight <- E(net)$raw_count

#Apply Disparity Filter [Serrano, Boguna & Vespignani (2009) PNAS 106(16)]
bb <- disparity_filter(net = net)

rm(net)

#Analyse the applied filter
rem <- analyse_filter(bb)

#Generate a filtered network with automated disparity threshold (best trade off between fraction of remaining nodes and links)
filtered.data <- recommend_threshold(bb, return.network = T)

#Plot results
#Fraction of remaining links vs nodes
plots[[1]] <- ggplot(rem, aes(L, N, colour = filter)) + geom_line(size = 1.5) + scale_x_reverse() + 
  annotate("point", x = filtered.data$threshold$L, y = filtered.data$threshold$N, colour = "black", size = 3) + 
  labs(x = expression(L[bb]/L[tot]), y = expression(N[bb]/N[tot])) + theme_bw() + 
  theme(legend.title = element_blank(), legend.background = element_blank(), legend.justification=c(0,0), legend.position=c(0,0))

#Fraction of remaining total weight vs nodes
plots[[2]] <- ggplot(rem, aes(W, N, colour = filter)) + geom_line(size = 1.5) + scale_x_reverse() + 
  labs(x = expression(W[bb]/W[tot]), y = expression(N[bb]/N[tot])) + theme_bw() + 
  theme(legend.title = element_blank(), legend.background = element_blank(), legend.justification=c(0,0), legend.position=c(0,0))

#Resulting node degree distribution
plots[[3]] <- plot_degree_distr(network = filtered.data$filtered.net)

##Fraction of remaining links vs fraction of nodes in LCC
plots[[4]] <- ggplot(rem, aes(L, LCC, colour = filter)) + geom_point(size = 1.5) + scale_x_reverse() + 
  labs(x = expression(L[bb]/L[tot]), y = "Fraction of nodes in LCC") + theme_bw() + 
  theme(legend.title = element_blank(), legend.background = element_blank(), legend.justification=c(0,0), legend.position=c(0,0))

plots.disparity <- plot_grid(plotlist = plots, nrow = 2, ncol = 2, labels = letters[1:4])
save_plot("../results/disparity_filter.pdf", plots.disparity, nrow = 2, ncol = 2, base_aspect_ratio = 1.3)

save(bb, rem, file = "../results/disparity_net.RData")
