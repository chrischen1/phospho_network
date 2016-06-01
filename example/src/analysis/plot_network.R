output_path = '~/Documents/workspace/phospho_network/example/script_files/network_plot/'
script_path = '~/Documents/workspace/phospho_network/example/script_files/output'
network_file = '/home/cc400/Documents/workspace/temp_files/network.csv'
dist_cutoff  = 4
impo_cutoff  = 0.95
max_remot_num = 20
configs = c('enL_heat_results','enL_all_results','enL_heat_nop_results','svmL_heat_results','svmL_all_results','rfL_heat_results','rfL_all_results')

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
  selected_nodes <- beta_result$predictor[abs(apply(beta_result[,-(1:2)],1,median)) > max(0,quantile(abs(apply(beta_result[,-(1:2)],1,median)),impo_cutoff))]
  predictor_nodes <- unique(c(test_node_valid,intersect(colnames(nodes_d)[apply(nodes_d,2,max) <= dist_cutoff],selected_nodes)))
  sp <- get.shortest.paths(network, from=test_node_valid, to=predictor_nodes,mode = 'all',output="vpath")
  sub_nodes <- unlist(sp$vpath,use.names = F)
  sub_g <- induced_subgraph(network, sub_nodes)
  node_imp <- apply(beta_result[,-(1:2)],1,median)
  names(node_imp) <- beta_result$predictor
  node_imp_nonz <- node_imp[node_imp > 0]
  node_imp_valid <- node_imp_nonz[node_imp_nonz > quantile(node_imp_nonz,0.95)]
  node_imp_add <- node_imp_valid[!names(node_imp_valid) %in% V(sub_g)$name]
  node_imp_add <- sort(node_imp_add,decreasing = T)[1:min(length(node_imp_add),max_remot_num)]
  node_imp_add <- node_imp_add[!is.na(node_imp_add)]
  sub_g <- sub_g + vertices(names(node_imp_add))+vertices('remote_nodes')
  node_imp_g <- abs(node_imp[names(node_imp)[names(node_imp) %in% V(sub_g)$name]])/sum(abs(node_imp[names(node_imp)[names(node_imp) %in% V(sub_g)$name]]))
  sub_g_size <- rep(1,length(V(sub_g)$name))
  names(sub_g_size) <- V(sub_g)$name
  sub_g_size[names(node_imp_g[node_imp_g>0])] <- node_imp_g[node_imp_g>0]/min(node_imp_g[node_imp_g>0]) *2
  for(rnode in names(node_imp_add)){
    sub_g <- add.edges(sub_g,c(rnode,'remote_nodes'))
  }
  for(snode in test_node_valid){
    sub_g <- add.edges(sub_g,c('remote_nodes',snode))
  }
  V(sub_g)$size <- sub_g_size/max(sub_g_size)*10+4
  V(sub_g)$color[V(sub_g)$name %in% names(node_imp_add)] <- 'orange'
  V(sub_g)$color[V(sub_g)$name %in% predictor_nodes] <- 'red'
  V(sub_g)$color[V(sub_g)$name %in% test_node_valid] <- 'green'
  V(sub_g)$label.cex = 0.8
  plot(sub_g,edge.arrow.size=0.4)
}
plot_config <- function(script_path,config,output_path,network){
  beta_file_list   <- grep('beta',list.files(paste(script_path,config,sep = '/')),value = T)
  beta_table <- NULL
  rp <- paste(script_path,config,sep = '/')
  dir.create(paste(output_path,config,sep = ''))
  for(beta_file in beta_file_list){
    if(is.null(beta)){
      beta_table <- read.csv(paste(rp,beta_file,sep = '/'),as.is = T)
    }else{
      beta_table <- rbind(beta_table,read.csv(paste(rp,beta_file,sep = '/'),as.is = T))
    }
  }
  test_names <- unique(beta_table$gene_site)
  for(test_name in test_names){
    test_node_valid <- target_select(test_name,network)
    beta_result <- beta_table[beta_table$gene_site == test_name,]
    if(length(test_node_valid) > 0){
      print(test_name)
      png(filename = paste(output_path,config,'/',test_node_valid,'.png',sep = ''),width = 1280,height = 960)
      plot_pred(test_node_valid,beta_result,network)
      dev.off()
    }
  }
}

network <- graph_from_adjacency_matrix(as.matrix(read.csv(network_file,row.names = 1)),mode = 'directed')
for(config in configs){
  plot_config(script_path,config,output_path,network)
}


