require(igraph)
require(MASS)


el2adj <- function(el,d = F){
  adj <- as.matrix(get.adjacency(graph_from_edgelist(as.matrix(el),directed = d)))
  return(adj)
}

heat_influence <- function(adj,b=0.6){
  w <- t(t(adj)/colSums(adj))
  influence_matrix <- ginv(diag(nrow(adj))-(1-b)*w)*b
  colnames(influence_matrix) <- colnames(adj)
  rownames(influence_matrix) <- rownames(adj)
  return(t(influence_matrix))
}

heat_difussion_analysis <- function(adj,influence_matrix,y_ind){
  net_im=graph.adjacency(t(influence_matrix),mode="directed",weighted=TRUE,diag=F) #"heat" should flow from downstream of the pathway to upstream
  g <- graph_from_adjacency_matrix(adj)
  plot(g,vertex.label.cex = 1,edge.arrow.size=0,vertex.size = 6.3)
  # plot influence matrix
  E <- t(apply(get.edgelist(net_im),1,sort))
  E(net_im)$curved <- 0
  E(net_im)[duplicated(E) | duplicated(E,fromLast =TRUE)]$curved <- 0.4
  w_vec <- (log(E(net_im)$weight+0.1)-min(log(E(net_im)$weight+0.1)))*3
  plot(net_im,edge.width=w_vec,vertex.label.cex = 0.4,edge.arrow.size = max(w_vec)/10)
  # print(influence_matrix)
  #heat diffusion
  heat_matrix <- matrix(0,ncol = ncol(influence_matrix),nrow = nrow(influence_matrix))
  colnames(heat_matrix) <- colnames(influence_matrix)
  rownames(heat_matrix) <- rownames(influence_matrix)
  heat_matrix[y_ind,y_ind] <- 1
  heat_diffusion <- influence_matrix %*% heat_matrix
  net_hd=graph.adjacency(t(heat_diffusion),mode="directed",weighted=TRUE,diag=F)
  # plot(net_hd,edge.width=E(net_hd)$weight*7/max(E(net_hd)$weight),edge.label=round(E(net_hd)$weight,3),vertex.label.cex = 0.7)
  # print(heat_diffusion)
  plot_heat_flow(network_adj = adj,influence_matrix = influence_matrix,predictors = heat_diffusion[heat_diffusion[,y_ind]>0.01*max(heat_diffusion[,y_ind]),y_ind])
}

add_gfeatures <- function(edge_list){
  node_index <- grep('psite',unique(c(edge_list[,1],edge_list[,2])),value = T,invert = T)
  cnvs <- paste(node_index,'CNV',sep = '_')
  muts <- paste(node_index,'mutation',sep = '_')
  mrna <- paste(node_index,'RNA',sep = '_')
  new_el <- rbind(edge_list,cbind(c(muts,mrna),rep(node_index,2)),cbind(cnvs,mrna))
  return(new_el)
}

predictor_construct <- function(test_node,test_gene,influence_matrix,b = 0.6, k = 0.05){
  #create heat flow by heat influence and test_node
  heat_input <- matrix(0,nrow(influence_matrix),nrow(influence_matrix))
  colnames(heat_input) <- colnames(influence_matrix)
  rownames(heat_input) <- colnames(influence_matrix)
  heat_input[test_node,test_node] <- 1
  heat_flow <- influence_matrix %*% heat_input
  predictors <- heat_flow[,test_node][heat_flow[,test_node] >= k*heat_flow[paste(test_gene,'CNV',sep = '_'),test_node]]
  return(predictors)
}

# plot predictors with heat flow in orginal pathway structure
plot_heat_flow <- function(network_adj,influence_matrix,predictors,target_name){
  #create sub heat influence matrix by predictors
  influence_matrix2 <- influence_matrix[names(predictors),names(predictors)]
  influence_matrix2[influence_matrix2<0] <- 0
  influence_matrix2 <- t(influence_matrix2)
  influence_matrix2[t(network_adj[names(predictors),names(predictors)])==0] <- 0
  #normalize heat value and sub heat influence matrix
  predictors <- predictors/sum(predictors)
  percent_heat_receive <- influence_matrix2/apply(influence_matrix2,1,sum)
  percent_heat_receive[is.nan(percent_heat_receive)] <- 0
  # each node recevie heat by its heat from predictors times the ratio of heat receive in sub heat influence matrix
  g_sub_adj <- percent_heat_receive * predictors
  for (i in 1:(nrow(g_sub_adj)-1)){
    for(j in (i+1):nrow(g_sub_adj)){
      if(g_sub_adj[i,j]<g_sub_adj[j,i]){
        g_sub_adj[i,j] <- 0
      }else{
        g_sub_adj[j,i] <- 0
      }
    }
  }
  g_sub <- graph_from_adjacency_matrix(g_sub_adj,weighted = T)
  w_vec <- log(E(g_sub)$weight)-min(log(E(g_sub)$weight))+0.1
  V(g_sub)$color <- 'white'
  V(g_sub)$color[grep('_',V(g_sub)$name)] <- 'orange'
  V(g_sub)$color[which(V(g_sub)$name == target_name)] <- 'red'
  plot(g_sub,edge.width=w_vec*1.1,edge.arrow.size = max(w_vec)/66,vertex.label.cex = 1,vertex.size = 7.3)
}

#fig2a
edge_list <- NULL
for (i in 1:4){
  edge_list <- rbind(edge_list,cbind(paste('Gene',i,sep = ''),paste('Gene',i+1,sep = '')))
}
edge_list <- rbind(edge_list,cbind('Gene3',paste('Gene3','psite',sep = '_')))
edge_list <- rbind(edge_list,cbind('Gene2',paste('Gene3','psite',sep = '_')))
plot(graph_from_edgelist(edge_list))

#fig2b
new_network <- add_gfeatures(edge_list)
target_name <- 'Gene3_psite'
adj <- el2adj(new_network,d = F)
influence_matrix <- heat_influence(t(adj))
net_im=graph.adjacency(t(influence_matrix),mode="directed",weighted=TRUE,diag=F) #"heat" should flow from downstream of the pathway to upstream
g <- graph_from_adjacency_matrix(adj,mode = 'undirected')
y_ind <- which(V(g)$name==target_name)
V(g)$color <- 'white'
V(g)$color[grep('_',V(g)$name)] <- 'orange'
V(g)$color[y_ind] <- 'red'
E(g)$curved <- 0.4
plot(g,vertex.label.cex = 1,edge.arrow.size=0,vertex.size = 8.3)

#fig2c
E <- t(apply(get.edgelist(net_im),1,sort))
E(net_im)$curved <- 0
E(net_im)[duplicated(E) | duplicated(E,fromLast =TRUE)]$curved <- 0.4
w_vec <- (log(E(net_im)$weight+0.1)-min(log(E(net_im)$weight+0.1)))*3
# plot(net_im,edge.width=w_vec,vertex.label.cex = 0.4,edge.arrow.size = max(w_vec)/10)
heat_matrix <- matrix(0,ncol = ncol(influence_matrix),nrow = nrow(influence_matrix))
colnames(heat_matrix) <- colnames(influence_matrix)
rownames(heat_matrix) <- rownames(influence_matrix)
heat_matrix[y_ind,y_ind] <- 1
heat_diffusion <- influence_matrix %*% heat_matrix
net_hd=graph.adjacency(t(heat_diffusion),mode="directed",weighted=TRUE,diag=F)
plot_heat_flow(network_adj = adj,influence_matrix,predictors = heat_diffusion[,y_ind],target_name)

