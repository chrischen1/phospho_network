# This Rscript would take a processed measured data file as input, and do the following steps:
# 1.Import the protein data and network file(as edge list) from csv files
# 2.Create a heat influence matrix by HOTNET
# input files:
# 1. rna/cnv_filename and missense_table_file/missense_table_file are processed data for predictors
# 2. interaction_gene_file: info for gene interactions or gene interactions to build the network
# 3. interaction_site_file: addtional info for protein/site interactions or gene interactions
# output files:
# 1.heat_outfile: heat influence matrix
# 2.network_outfile network for prediction

library(methods)
library(igraph)
#library(MASS)

# Set inputs:
rna_filename          = '~/Documents/workspace/phospho_network/example/script_files/rna_processed.csv'
cnv_filename          = '~/Documents/workspace/phospho_network/example/script_files/cnv_processed.csv'
mutation_matrix_file  = '~/Documents/workspace/phospho_network/example/script_files/mutation_matrix_misc.csv'
missense_table_file   = '~/Documents/workspace/phospho_network/example/script_files/mutation_missense.csv'
interaction_site_file = '~/Documents/workspace/phospho_network/example/script_files/network/interaction_table_site_all.csv'
interaction_gene_file = '~/Documents/workspace/phospho_network/example/script_files/network/interaction_table_gene_all.csv'

# Set parameters
heat_outfile   = '~/Documents/workspace/temp_files/heat_influence_ud.csv' 
network_outfile   = '~/Documents/workspace/temp_files/network_ud.csv'

b = 0.6  # parameter to define the effciency of heat diffusion
d = F    # the heat flow on graph is directed?
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

# Main
# import dataset
rna_data  <- read.csv(rna_filename,as.is = T,row.names = 1)
cnv_data  <- read.csv(cnv_filename,as.is = T,row.names = 1)
mutation_misc  <- read.csv(mutation_matrix_file,as.is = T,row.names = 1)
missense_table <- read.csv(missense_table_file,as.is = T,row.names = 1) 

# construct network
site_interactions <- data.frame(read.csv(interaction_site_file,as.is = T))
gene_interactions <- data.frame(read.csv(interaction_gene_file,as.is = T))

# keep only the intersect of RNA measured and interactions
gene_network   <- unique(c(gene_interactions$GeneA,gene_interactions$geneB)) #use this on HPC
# gene_network   <- unique(c(site_interactions$GeneA,site_interactions$geneB)) #use this for test

gene_measured  <- unique(c(rownames(rna_data),rownames(cnv_data),rownames(mutation_misc),missense_table$gene))
gene_intersect <- intersect(gene_network,gene_measured)

gene_mut <- intersect(rownames(mutation_misc),gene_intersect)
gene_cnv <- intersect(rownames(cnv_data),gene_intersect)
gene_rna <- unique(c(intersect(rownames(rna_data),gene_intersect),gene_cnv))
gene_mis <- intersect(missense_table$gene,gene_intersect)
gene_interactions_intersect <- gene_interactions[gene_interactions$GeneA %in% gene_intersect & gene_interactions$GeneB %in% gene_intersect,]
site_interactions_intersect <- site_interactions[site_interactions$GeneA %in% gene_intersect & site_interactions$GeneB %in% gene_intersect,]
missense_intersect <- missense_table[missense_table$gene %in% gene_mis,]

edge_list <- matrix(NA,nrow = 0,ncol = 2)
colnames(edge_list) <- c('nodeA','nodeB')
edge_list <- rbind(edge_list,cbind(gene_interactions_intersect$GeneA,gene_interactions_intersect$GeneB))
edge_list <- rbind(edge_list,cbind(site_interactions_intersect$GeneB,paste(site_interactions_intersect$GeneB,site_interactions_intersect$Position,sep = '_')))
edge_list <- rbind(edge_list,cbind(site_interactions_intersect$GeneA,paste(site_interactions_intersect$GeneB,site_interactions_intersect$Position,sep = '_')))
edge_list <- rbind(edge_list,cbind(paste(gene_mis,'mutation_misc',sep = '_'),gene_mis))
edge_list <- rbind(edge_list,cbind(paste(gene_rna,'RNA',sep = '_'),gene_rna))
edge_list <- rbind(edge_list,cbind(paste(gene_cnv,'CNV',sep = '_'),paste(gene_cnv,'RNA',sep = '_')))
edge_list <- rbind(edge_list,cbind(missense_intersect$missense_names,missense_intersect$gene))
edge_list <- as.matrix(unique(apply(edge_list,2,as.character)))


network_data <- graph_from_edgelist(edge_list,directed = d)
network_adj  <- as.matrix(get.adjacency(network_data))

#create heat influence matrix by adjacency matrix
hi_matrix <- heat_influence(t(network_adj),b = b)
write.csv(hi_matrix,heat_outfile)
write.csv(network_adj,network_outfile)
