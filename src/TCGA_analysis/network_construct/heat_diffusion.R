# This Rscript would take a processed measured data file as input, and do the following steps:
# 1.Import the protein data and network file(as edge list) from csv files
# 2.Create a heat influence matrix by HOTNET
# input files:
# 1. measure_filename are processed quantitive data for phos protein
# 2. interaction_gene_file: info for gene interactions or gene interactions to build the network
# 3. interaction_site_file: addtional info for protein/site interactions or gene interactions

# output files:
# 1.heat_outfile: heat influence matrix
# 2.network_outfile network for prediction

library(methods)
library(igraph)
library(MASS)

# Set inputs:
rna_filename          = '~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/BRCA.exp.547.med.txt'
interaction_site_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/interaction_table_site_all.csv'
interaction_gene_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/interaction_table_gene_all.csv'

# Special for RPPA analysis
measure_filename      = '~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/rppaData-403Samp-171Ab-Trimmed.txt'
rppa_mapping_file1    = '~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/BRCA.RPPA.Level_3/mdanderson.org_BRCA.MDA_RPPA_Core.mage-tab.1.0.0/mdanderson.org_BRCA.MDA_RPPA_Core.antibody_annotation.txt'
rppa_mapping_file2    = '~/Documents/workspace/phospho_network/RAWDATA/TCPA/TCPA_ANTIBODY_MAPPING.txt'

# Set parameters
heat_outfile          = '~/Documents/workspace/phospho_network/script_files/TCGA_analysis/heat_influence.csv'
network_outfile       = '~/Documents/workspace/phospho_network/script_files/TCGA_analysis/network.csv'


b = 0.6  # a parameter to define the effciency of heat diffusion


# Functions:
#add additional edges based on gene-gene interactions and gene-site interactions
add_edges <- function(gene_interactions,site_interactions){
  gene_list <- unique(c(gene_interactions$GeneA,gene_interactions$GeneB))
  edge_list <- matrix(NA,nrow = 0,ncol = 2)
  colnames(edge_list) <- c('nodeA','nodeB')
  edge_list <- rbind(edge_list,cbind(gene_interactions$GeneA,gene_interactions$GeneB))
  edge_list <- rbind(edge_list,cbind(site_interactions$GeneB,paste(site_interactions$GeneB,site_interactions$Position,sep = '_')))
  edge_list <- rbind(edge_list,cbind(site_interactions$GeneA,paste(site_interactions$GeneB,site_interactions$Position,sep = '_')))
  edge_list <- rbind(edge_list,cbind(paste(gene_list,'mutation',sep = '_'),gene_list))
  edge_list <- rbind(edge_list,cbind(paste(gene_list,'RNA',sep = '_'),gene_list))
  edge_list <- rbind(edge_list,cbind(paste(gene_list,'CNV',sep = '_'),paste(gene_list,'RNA',sep = '_')))
  edge_list_final <- as.matrix(unique(apply(edge_list,2,as.character)))
  return(edge_list_final)
}

#hotnet2 heat diffusion
heat_influence <- function(adj,b=0.6){
  w <- t(t(adj)/colSums(adj))
  w[adj==0] <- 0
  influence_matrix <- ginv(diag(nrow(adj))-(1-b)*w)*b
  colnames(influence_matrix) <- colnames(adj)
  rownames(influence_matrix) <- rownames(adj)
  return(t(influence_matrix))
}

# Main
# import dataset
rna_data_raw  <- read.table(rna_filename,as.is = T,row.names = 1,fill = NA)
rna_data <- rna_data_raw[-1,]
colnames(rna_data) <- gsub('-','.',rna_data_raw[1,-1],fixed = T)
colnames(rna_data) <- gsub(pattern = '(\\w+\\.\\w+\\.\\w+\\.\\d+).+',replacement = '\\1',x = colnames(rna_data))

# construct network
site_interactions <- data.frame(read.csv(interaction_site_file,as.is = T))
gene_interactions <- data.frame(read.csv(interaction_gene_file,as.is = T))

# keep only the intersect of RNA measured and interactions
# gene_list <- unique(rownames(rna_data)) # would use this on HPC
gene_list <- intersect(rownames(rna_data),unique(c(site_interactions$GeneA,site_interactions$ProteinB)))
genes_intersect <- intersect(gene_list,unique(c(gene_interactions$GeneA,gene_interactions$geneB)))
gene_interactions_intersect <- gene_interactions[gene_interactions$GeneA %in% genes_intersect & gene_interactions$GeneB %in% genes_intersect,]
site_interactions_intersect <- site_interactions[site_interactions$GeneA %in% genes_intersect & site_interactions$GeneB %in% genes_intersect,]

edge_list <- add_edges(gene_interactions_intersect,site_interactions_intersect)
network_data <- graph_from_edgelist(edge_list,directed = T)
network_adj  <- as.matrix(get.adjacency(network_data))

#create heat influence matrix by adjacency matrix
hi_matrix <- heat_influence(t(network_adj),b = b)
write.csv(hi_matrix,heat_outfile)
write.csv(network_adj,network_outfile)
