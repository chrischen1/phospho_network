require(igraph)
edge_list = matrix(c('B','A_psite1','B','A_psite2','C','A_psite2'),ncol = 2,byrow = T)
g <- graph_from_edgelist(edge_list)
V(g)$color <- 'orange'
V(g)$color[V(g)$name %in% c('A_psite1','A_psite2')] <- 'red'
plot(g)
