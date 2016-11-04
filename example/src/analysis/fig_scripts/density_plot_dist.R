library(igraph)
network_pred_el = read.csv('~/Documents/workspace/phospho_network/example/script_files/network/edge_list_ud.csv',as.is = T)
network_test_el  = read.csv('~/Documents/workspace/phospho_network/example/script_files/network/pw_table_ppi.csv',as.is = T)
source_node = c('GSK3B','GSK3A')
dist_cutoff = 5
sig_nodes = c('RAB5C','MMP26','PXDNL')

network_pred <- graph_from_edgelist(as.matrix(network_pred_el))
network_test <- graph_from_edgelist(as.matrix(network_test_el))
dist_nodes <- apply(distances(network_pred,source_node),2,min)
test_nodes <- names(dist_nodes)[dist_nodes>dist_cutoff-1]
test_genes <- intersect(test_nodes,V(network_test)$name)

dist_nodes_large <- apply(distances(network_test,source_node,test_genes),2,min)
dist_nodes_large[is.infinite(dist_nodes_large)] <- max(dist_nodes_large[!is.infinite(dist_nodes_large)])
dist_nodes_sig <- apply(distances(network_test,source_node,sig_nodes),2,min)
