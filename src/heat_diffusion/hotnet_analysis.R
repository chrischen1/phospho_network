require(igraph)
require(MASS)

el2adj <- function(el,d = F){
  adj <- as.matrix(get.adjacency(graph_from_edgelist(as.matrix(el),directed = d)))
  return(adj)
}

heat_influence <- function(adj,b=0.6){
  w <- t(t(adj)/colSums(adj))
  w[adj==0] <- 0
  influence_matrix <- ginv(diag(nrow(adj))-(1-b)*w)*b
  colnames(influence_matrix) <- colnames(adj)
  rownames(influence_matrix) <- rownames(adj)
  return(influence_matrix)
}

heat_difussion_analysis <- function(edge_list,node_index,influence_matrix,y_ind = 1){
  node_names <- node_index$V2
  colnames(influence_matrix) <- node_names
  rownames(influence_matrix) <- node_names
  net_im=graph.adjacency(t(influence_matrix),mode="directed",weighted=TRUE,diag=F) #"heat" should flow from downstream of the pathway to upstream
  edge_list <- as.matrix(edge_list)
  g <- graph_from_edgelist(edge_list,directed = T)
  g_adj <- as.matrix(get.adjacency(g))
  colnames(g_adj) <- node_names
  rownames(g_adj) <- node_names
  g <- graph.adjacency(g_adj,mode="directed")
  plot(g)
  
  # plot influence matrix
  E <- t(apply(get.edgelist(net_im),1,sort))
  E(net_im)$curved <- 0
  E(net_im)[duplicated(E) | duplicated(E,fromLast =TRUE)]$curved <- 0.2
  plot(net_im,edge.width=E(net_im)$weight*10/max(influence_matrix),edge.arrow.size = max(influence_matrix)*0.8)
  print(influence_matrix)
  #heat diffusion
  heat_matrix <- matrix(0,ncol = ncol(influence_matrix),nrow = nrow(influence_matrix))
  colnames(heat_matrix) <- colnames(influence_matrix)
  rownames(heat_matrix) <- rownames(influence_matrix)
  heat_ind <- which(rownames(heat_matrix) == node_index$V2[y_ind])
  heat_matrix[heat_ind,heat_ind] <- 1
  heat_diffusion <- t(influence_matrix %*% heat_matrix)
  net_hd=graph.adjacency(heat_diffusion,mode="directed",weighted=TRUE,diag=F)
  plot(net_hd,edge.width=E(net_hd)$weight*7/max(E(net_hd)$weight),edge.label=round(E(net_hd)$weight,3),edge.arrow.size = max(heat_diffusion)*0.8)
  print(heat_diffusion)
}

add_gfeatures <- function(edge_list,node_index){
  cnvs <- paste(node_index[,2],'cnv',sep = '_')
  muts <- paste(node_index[,2],'mut',sep = '_')
  mrna <- paste(node_index[,2],'mrna',sep = '_')
  newf <- cbind((max(node_index[,1])+1):(max(node_index[,1])+length(c(cnvs,muts,mrna))),c(cnvs,muts,mrna))
  new_ni <- rbind(node_index,newf)
  newi <- cbind(newf[,1],node_index[,1])
  new_el <- rbind(edge_list,apply(newi,2,as.numeric))
  return(list('node_index'=new_ni,'edge_list'= new_el))
}

#main body
# a toy model:
edge_list <- read.table('~/Downloads/hotnet2/toy_model/example_edgelist.txt',header = F,sep = ' ',as.is = T)
node_index <- read.table('~/Downloads/hotnet2/toy_model/example_gene_index.txt',header = F,sep = ' ',as.is = T)
adj <- el2adj(edge_list,d = T)
influence_matrix <- heat_influence(adj)
heat_difussion_analysis(edge_list,node_index,influence_matrix,y_ind = 3)

#Real model: MAPK pathway
mapk_adj <- read.csv('~/Documents/workspace/phospho_network/RAWDATA/mapk_info/mapk_network_kegg.csv',row.names = 1,as.is = T)
mapk_adj_num <- as.matrix(mapk_adj)
colnames(mapk_adj_num) <- 1:nrow(mapk_adj_num)
rownames(mapk_adj_num) <- 1:nrow(mapk_adj_num)
node_index_mapk <- cbind(colnames(mapk_adj_num),colnames(mapk_adj))
node_index_mapk[,1] <- as.numeric(node_index_mapk[,1])
node_index_mapk <- as.data.frame(node_index_mapk)
edge_list_mapk <- get.edgelist(graph.adjacency(mapk_adj_num,mode="directed"))
edge_list_mapk <- apply(edge_list_mapk,2,as.numeric)
influence_matrix_mapk <- heat_influence(el2adj(edge_list_mapk,d = T))
heat_difussion_analysis(edge_list_mapk,node_index_mapk,influence_matrix_mapk,y_ind = which(node_index_mapk$V2 == 'RAF1'))

#add genetic features for simple model:
new_network <- add_gfeatures(edge_list,node_index)
edge_list2 <- new_network$edge_list
node_index2 <- new_network$node_index
adj2 <- el2adj(edge_list2,d = T)
influence_matrix2 <- heat_influence(adj2)
heat_difussion_analysis(edge_list2,node_index2,influence_matrix2,y_ind = 3)
