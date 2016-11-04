ggi_file = '~/Documents/workspace/phospho_network/script_files/network/interaction_table_gene_all.csv'
gsi_file = '~/Documents/workspace/phospho_network/script_files/network/interaction_table_site_all.csv'

ms_p_file = '~/Documents/workspace/phospho_network/script_files/ms_data_processed/msdata_processed.csv'
ms_total_file = '~/Documents/workspace/phospho_network/script_files/ms_data_processed/total_protein_processed.csv'
mut_file = '~/Documents/workspace/phospho_network/script_files/ms_data_processed/mutation_matrix.csv'

output_path = '~/Documents/workspace/temp_files'

library(igraph)
#hotnet2 heat diffusion
heat_influence <- function(adj,b=0.6){
  w <- t(t(adj)/colSums(adj))
  w[is.nan(w)] <- 0
  w[adj==0] <- 0
  influence_matrix <- solve(diag(nrow(adj))-(1-b)*w)*b
  colnames(influence_matrix) <- colnames(adj)
  rownames(influence_matrix) <- rownames(adj)
  return(t(influence_matrix))
}

ggi <- read.csv(ggi_file,as.is = T)
gsi <- read.csv(gsi_file,as.is = T)
ms_p <- read.csv(ms_p_file,as.is = T)
ms_tot <- read.csv(ms_total_file,as.is = T)
mut_matrix <- read.csv(mut_file,row.names = 1)

mut_genes <- rownames(mut_matrix)[apply(mut_matrix,1,sum)>1]
network_genes <- intersect(unlist(ggi),c(ms_p$gene_symbol,ms_tot$Gene.Symbol,mut_genes))

gsi2 <- gsi[gsi$GeneA %in% network_genes & gsi$GeneB %in% network_genes,]
ms_p2 <- ms_p[ms_p$gene_symbol %in% network_genes,]
ms_tot2 <- ms_tot[ms_tot$Gene.Symbol %in% network_genes,]

ggi_edge <- ggi[ggi$GeneA %in% network_genes & ggi$GeneB %in% network_genes,]
gsi_edge <- cbind(gsi2$GeneA,paste(gsi2$GeneB,gsi2$Position,sep = '_'))
colnames(gsi_edge) <- colnames(ggi_edge)
gsi_edge2 <- cbind(gsi2$GeneB,paste(gsi2$GeneB,gsi2$Position,sep = '_'))
colnames(gsi_edge2) <- colnames(ggi_edge)
mut_edge <- cbind(paste(mut_genes,'mutation_misc',sep = '_'),mut_genes)
colnames(mut_edge) <- colnames(ggi_edge)
tot_edge <- (cbind(paste(ms_tot2$Gene.Symbol,'total',sep = '_'),ms_tot2$Gene.Symbol))
colnames(tot_edge) <- colnames(ggi_edge)
network_edgelist <- rbind(ggi_edge,gsi_edge,gsi_edge2,mut_edge,tot_edge)
network_edgelist <- as.matrix(unique(network_edgelist[network_edgelist$GeneA != '' & network_edgelist$GeneB != '',]))
network <- graph_from_edgelist(network_edgelist)

adj <- get.adjacency(network,type = 'both',sparse = F)
hi_matrix <- heat_influence(t(adj))


write.csv(adj,paste(output_path,'network_ms.csv',sep = '/'))
write.csv(hi_matrix,paste(output_path,'heat_matrix_ms.csv',sep = '/'))

