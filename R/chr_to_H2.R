
library(NetHypGeom)

load("results/disp_analysis.RData")

# Disparity-filter net
df_net <- delete_edges(net, E(net)[edge.pvals > 
                                     analysis$threshold[analysis$recommended]])
df_net <- delete_vertices(df_net, which(degree(df_net) < 1))

# Chromosome of interest
chr <- "chr1"

chr_net <- induced_subgraph(df_net, vids = which(V(df_net)$chr == chr))
chr_net <- delete_vertices(chr_net, which(degree(chr_net) < 1))

comp <- components(chr_net)
chr_net <- induced_subgraph(chr_net, vids = which(comp$membership == 
                                                    which.max(comp$csize)))

# Determine net properties, including Temperature
N <- vcount(chr_net)
avg_k <- mean(degree(chr_net))
gma = 3
cc = transitivity(chr_net, "average")

clust.at.zero <- numeric(10)

for(i in 1:length(clust.at.zero)){
  ps <- ps_model(N = N, avg.k = avg_k, gma = gma, Temp = 0)
  clust.at.zero[i] <- transitivity(ps$network, "average")
}

slope <- (0 - 1)/mean(clust.at.zero)
Temp <- slope*cc + 1

# Map net to hyperbolic space using LaBNE+HM
coords <- labne_hm(chr_net, gma = gma, Temp = Temp, w = 0)

# Visualise mapping
plot_hyperbolic_net(network = chr_net, nodes = coords$polar, 
                    node.colour = coords$polar$theta)
