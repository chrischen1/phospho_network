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
  return(t(influence_matrix))
}

heat_difussion_analysis <- function(adj,influence_matrix,y_ind){
  net_im=graph.adjacency(t(influence_matrix),mode="directed",weighted=TRUE,diag=F) #"heat" should flow from downstream of the pathway to upstream
  g <- graph_from_adjacency_matrix(adj)
  # plot(g,vertex.label.cex = 1,edge.arrow.size=0.2,vertex.size = 6.3)
  # plot influence matrix
  # E <- t(apply(get.edgelist(net_im),1,sort))
  # E(net_im)$curved <- 0
  # E(net_im)[duplicated(E) | duplicated(E,fromLast =TRUE)]$curved <- 0.4
  # w_vec <- (log(E(net_im)$weight+0.1)-min(log(E(net_im)$weight+0.1)))*3
  # plot(net_im,edge.width=w_vec,vertex.label.cex = 0.4,edge.arrow.size = max(w_vec)/10)
  # print(influence_matrix)
  #heat diffusion
  heat_matrix <- matrix(0,ncol = ncol(influence_matrix),nrow = nrow(influence_matrix))
  colnames(heat_matrix) <- colnames(influence_matrix)
  rownames(heat_matrix) <- rownames(influence_matrix)
  heat_ind <- which(rownames(heat_matrix) == y_ind)
  heat_matrix[heat_ind,heat_ind] <- 1
  heat_diffusion <- influence_matrix %*% heat_matrix
  net_hd=graph.adjacency(t(heat_diffusion),mode="directed",weighted=TRUE,diag=F)
  # plot(net_hd,edge.width=E(net_hd)$weight*7/max(E(net_hd)$weight),edge.label=round(E(net_hd)$weight,3),vertex.label.cex = 0.7)
  predictors = heat_diffusion[heat_diffusion[,heat_ind]>0.01*max(heat_diffusion[,heat_ind]),heat_ind]
  plot_heat_flow(network_adj = adj,influence_matrix,predictors,y_ind)
}

add_gfeatures <- function(edge_list){
  node_index <- unique(c(edge_list[,1],edge_list[,2]))
  cnvs <- paste(node_index,'cnv',sep = '_')
  muts <- paste(node_index,'mutation',sep = '_')
  mrna <- paste(node_index,'mrna',sep = '_')
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
plot_heat_flow <- function(network_adj,influence_matrix,predictors,y_ind = names((which.max(predictors)))){
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
  g_sub <- graph_from_adjacency_matrix(g_sub_adj,weighted = T)
  w_vec <- log(E(g_sub)$weight,2)-min(log(E(g_sub)$weight),2)+0.1
  V(g_sub)$color <- 'white'
  V(g_sub)$color[V(g_sub)$name %in% y_ind] <- 'red'
  V(g_sub)$name <- paste('         ',V(g_sub)$name,format(round(predictors[V(g_sub)$name],3),nsmall = 3))
  plot(g_sub,edge.width=w_vec,edge.arrow.size = 0,vertex.size = 10.2,vertex.label.cex = 1.6)
}
par(mfrow=c(1,3))
el_a <- matrix(c('A','C','B','C','C','S','E','S','D','E','F','E','H','S','G','H','I','H'),ncol = 2,byrow = T)
adj <- el2adj(el_a,d = F)
influence_matrix <- heat_influence(t(adj))
heat_difussion_analysis(adj,influence_matrix,y_ind = 'S')

el_c <- matrix(c('A','S1','B','S1','A','S2'),ncol = 2,byrow = T)
adj <- el2adj(el_c,d = F)
influence_matrix <- heat_influence(t(adj))
heat_difussion_analysis(adj,influence_matrix,y_ind = c('S1','S2'))

el_b <- matrix(c('S','A','S','B','S','C','S','D','A','E','B','F','C','G','D','G','H','G','I','F','E','G'),ncol = 2,byrow = T)
adj <- el2adj(el_b,d = F)
influence_matrix <- heat_influence(t(adj))
heat_difussion_analysis(adj,influence_matrix,y_ind = 'S')


