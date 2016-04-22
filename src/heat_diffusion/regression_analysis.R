# regression analysis based on hotnet2 hheat diffusion and gelnet package
# input: 
# 1.adj matrix for PPI
# 2.rna_expression, sample_mutation, total_protein
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

add_gfeatures <- function(edge_list){
  node_index <- unique(as.character(edge_list))
  cnvs <- paste(node_index,'cnv',sep = '_')
  muts <- paste(node_index,'mut',sep = '_')
  mrna <- paste(node_index,'mrna',sep = '_')
}


# main
mrna <- read.table('~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/BRCA.exp.547.med.txt',fill = NA)

  
  
  
  
