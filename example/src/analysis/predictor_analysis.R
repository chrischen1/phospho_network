args <- commandArgs(trailingOnly = TRUE)
test_name <- gsub('+',';',args[1],fixed = T)
config_filename <- args[2]
source(config_filename)

library(igraph,quietly = T)

pred_analysis <- function(test_name,network_adj){
  nodes <- colnames(network_adj)[colnames(network_adj) != test_name]
  predictors <- nodes[grep('_[A-Za-z]+',nodes)]
  pred_table <- matrix(0,ncol = 6,nrow = length(predictors))
  colnames(pred_table) <- c('target','predictor','type','dist_d','dist_ud','relation')
  pred_table[,1] <- test_name
  pred_table[,2] <- predictors
  pred_table[grep('CNV',predictors),3] <- 'CNV'
  pred_table[grep('RNA',predictors),3] <- 'RNA'
  pred_table[grep('RNA|CNV',predictors,invert = T),3] <- 'mutation'
  network <- graph_from_adjacency_matrix(network_adj,mode = 'directed')
  dist_d <- shortest.paths(network,v = test_name,to = predictors,algorithm = 'unweighted',mode = 'in')
  dist_ud <- shortest.paths(network,v = test_name,to = predictors,algorithm = 'unweighted')
  pred_table[,4] <- dist_d[,predictors]
  pred_table[,5] <- dist_ud[,predictors]
  pred_table[is.infinite(dist_ud[,predictors]),6] <- 'unconected'
  pred_table[!is.infinite(dist_d[,predictors]),6] <- 'upstream'
  pred_table[pred_table[,6]=='0',6] <- 'downstream'
  return(pred_table)
}

target_select <- function(mdata_names){
  mdata_p <- grep('\\w+_\\w+',mdata_names,value = T)
  mdata_plist <- strsplit(mdata_p,split = '_')
  mdata_pgenes <- strsplit(unlist(lapply(mdata_plist,function(x)x[1])),split = ';')
  mdata_psites <- strsplit(unlist(lapply(mdata_plist,function(x)x[2])),split = ';')
  mdata_pnodes <- c()
  for(i in 1:length(mdata_pgenes)){
    new_nodes=expand.grid(mdata_pgenes[[i]],mdata_psites[[i]])
    mdata_pnodes <- c(mdata_pnodes,paste(paste(new_nodes$Var1,new_nodes$Var2,sep = '_'),collapse = ';'))
  }
  return(cbind(mdata_p,mdata_pnodes))
}

print(test_name)
mdata <- read.csv(mdata_filename,row.names = 1)
network_adj   <- as.matrix(read.csv(network_file,row.names = 1))
node_all <- rownames(network_adj)
response_list <- target_select(rownames(mdata))

test_node <- strsplit(response_list[response_list[,1] == test_name,2],split = ';')[[1]]
test_node_valid <- test_node[test_node %in% node_all]

if(length(test_node_valid) == 0){
  test_node_valid <- unique(gsub('(\\w+)_\\w+','\\1',test_node))
  test_node_valid <- test_node_valid[test_node_valid %in% node_all]
}
if(length(test_node_valid) == 1){
  pred_table <- pred_analysis(test_node_valid,network_adj)
  pred_table[,1] <- test_name
  write.csv(pred_table,paste(result_path,paste(test_name,'pred_table.csv',sep = '_'),sep = '/'),row.names = F)
}else if(length(test_node_valid) > 1){
  pred_table_list <- list()
  for (i in 1:length(test_node_valid)){
    pred_table_list[[i]] <- pred_analysis(test_node_valid[i],network_adj)
  }
  pred_table <- pred_table_list[[1]]
  pred_table[,4] <- apply(do.call(cbind,lapply(pred_table_list,function(x)x[,4])),1,min)
  pred_table[,5] <- apply(do.call(cbind,lapply(pred_table_list,function(x)x[,5])),1,min)
  pred_table[,1] <- test_name
  write.csv(pred_table,paste(result_path,paste(test_name,'pred_table.csv',sep = '_'),sep = '/'),row.names = F)
}





