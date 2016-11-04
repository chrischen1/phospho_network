output_path = '~/Documents/workspace/phospho_network/example/script_files/network_plot/best_heat'
beta_file   = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_heat_best.csv'
network_file = '~/Documents/workspace/phospho_network/example/script_files/network/network_test.csv'
# network_file = '/home/cc400/Documents/workspace/temp_files/network.csv'
dist_cutoff  = 4
impo_cutoff  = 0.95
max_remot_num = 5
max_connect_num = 10
test_name <- 'EGFR_1068'

library(igraph)
target_select <- function(mdata_names,network){
  node_all <- V(network)$name
  mdata_p <- grep('\\w+_\\w+',mdata_names,value = T)
  mdata_plist <- strsplit(mdata_p,split = '_')
  mdata_pgenes <- strsplit(unlist(lapply(mdata_plist,function(x)x[1])),split = ';')
  mdata_psites <- strsplit(unlist(lapply(mdata_plist,function(x)x[2])),split = ';')
  mdata_pnodes <- c()
  for(i in 1:length(mdata_pgenes)){
    new_nodes=expand.grid(mdata_pgenes[[i]],mdata_psites[[i]])
    mdata_pnodes <- c(mdata_pnodes,paste(paste(new_nodes$Var1,new_nodes$Var2,sep = '_'),collapse = ';'))
  }
  response_list <- cbind(mdata_p,mdata_pnodes)
  test_node <- strsplit(response_list[response_list[,1] == mdata_names,2],split = ';')[[1]]
  test_node_valid <- test_node[test_node %in% node_all]
  if(length(test_node_valid) == 0){
    test_node_valid <- unique(gsub('(\\w+)_\\w+','\\1',test_node))
    test_node_valid <- test_node_valid[test_node_valid %in% node_all]
  }
  return(test_node_valid)
}
plot_pred <- function(test_node_valid,beta_result,network){
  nodes_d <- distances(network,v = test_node_valid)
  neighbor_nodes <- colnames(nodes_d)[apply(nodes_d,2,max) <= dist_cutoff]
  importance_nodes <- abs(apply(beta_result[,-(1:2)],1,median))
  names(importance_nodes) <- beta_result$predictor
  importance_nodes <- importance_nodes[importance_nodes>0]
  important_neighbor <- importance_nodes[names(importance_nodes) %in% neighbor_nodes]
  important_neighbor <- head(sort(important_neighbor,decreasing = T),max_connect_num)
  predictor_nodes <- names(important_neighbor)
  sub_nodes <- c()
  for(source_node in test_node_valid){
    sp <- get.shortest.paths(network, from=source_node, to=predictor_nodes,mode = 'all',output="vpath")
    sub_nodes <- unique(c(sub_nodes,unlist(sp$vpath,use.names = F)))
  }
  sub_g <- induced_subgraph(network, sub_nodes)
  node_imp_valid <- importance_nodes[importance_nodes > quantile(importance_nodes,impo_cutoff)]
  node_imp_remote <- node_imp_valid[!names(node_imp_valid) %in% V(sub_g)$name]
  node_imp_remote <- head(sort(node_imp_remote,decreasing = T),max_remot_num)
  node_imp_remote <- node_imp_remote[!is.na(node_imp_remote)]
  sub_g <- sub_g + vertices(names(node_imp_remote))+vertices('remote_nodes')
  node_imp <- c(important_neighbor,node_imp_remote)
  node_imp_g <- abs(node_imp[names(node_imp)[names(node_imp) %in% V(sub_g)$name]])/sum(abs(node_imp[names(node_imp)[names(node_imp) %in% V(sub_g)$name]]))
  sub_g_size <- rep(1,length(V(sub_g)$name))
  names(sub_g_size) <- V(sub_g)$name
  sub_g_size[names(node_imp_g[node_imp_g>0])] <- node_imp_g[node_imp_g>0]/min(node_imp_g[node_imp_g>0]) *2
  for(rnode in names(node_imp_remote)){
    sub_g <- add.edges(sub_g,c(rnode,'remote_nodes'))
  }
  for(snode in test_node_valid){
    sub_g <- add.edges(sub_g,c('remote_nodes',snode))
  }
  V(sub_g)$size <- sub_g_size/max(sub_g_size)*10+4
  V(sub_g)$color[V(sub_g)$name %in% names(node_imp_remote)] <- 'yellow'
  V(sub_g)$color[V(sub_g)$name %in% predictor_nodes] <- 'orange'
  V(sub_g)$color[V(sub_g)$name %in% test_node_valid] <- 'red'
  V(sub_g)$label.cex = 1
  write.table(as.matrix(as_adjacency_matrix(sub_g)),file = '~/Documents/workspace/phospho_network/egfr_sub_adj.tsv',sep = '\t')
  write.table(cbind(V(sub_g),V(sub_g)$size),file = '~/Documents/workspace/phospho_network/egfr_sub_v.tsv',sep = '\t')
  plot(sub_g,edge.arrow.size=0)
}

network <- graph_from_adjacency_matrix(as.matrix(read.csv(network_file,row.names = 1)),mode = 'undirected')
beta_table <- read.csv(beta_file,as.is = T)
test_node_valid <- target_select(test_name,network)
beta_result <- beta_table[beta_table$gene_site == test_name,]
pdf(file = '~/Documents/workspace/phospho_network/egfr_scatter.pdf')
plot_pred(test_node_valid,beta_result,network)
dev.off()