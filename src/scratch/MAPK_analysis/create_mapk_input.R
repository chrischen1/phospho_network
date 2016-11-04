# This script takes 4 input files:
# 1.interaction table for genes
# 2.prot interacitons table
# 3.site interacitons table
# 4.a list of genes in MAPK pathway which needs to predict
# the outputs are:
# 1.out_adj_m would be an adj table for genes based on literal
# 2.out_prot_table would be interaction table for porteins based on data from ks/signor
# 3.out_gene_table would be interaction table for porteins based on data from merged interaction table for genes from literal

library(igraph)
interaction_file_path      = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed'
mapk_info_path             = '~/Documents/workspace/phospho_network/RAWDATA/mapk_info'
mapk_result_path           = '~/Documents/workspace/phospho_network/script_files/mapk_analysis'

site_table_file  = 'interaction_table_site_all.csv'
prot_table_file  = 'interaction_table_prot_all.csv'
all_table_file   = 'network_all_gene.csv'
mapk_gene_file   = 'mapk_genes.txt'

out_adj_m        = 'mapk_network_literal.csv'
out_prot_table   = 'mapk_prot_table.csv'
out_gene_table   = 'mapk_gene_table.csv'

# Main
mapk_genes <- unlist(read.table(paste(mapk_info_path,mapk_gene_file,sep = '/'),as.is = T))

site_table <- read.csv(paste(interaction_file_path,site_table_file,sep = '/'),as.is = T)
prot_table <- read.csv(paste(interaction_file_path,prot_table_file,sep = '/'),as.is = T)
all_table  <- read.csv(paste(interaction_file_path,all_table_file,sep = '/'),as.is = T)

all_mapk     <- all_table[all_table$geneB %in% mapk_genes,]
all_mapk2    <- all_mapk[all_mapk$geneA != all_mapk$geneB,]
g_mapk       <- graph_from_edgelist(as.matrix(all_mapk2))
mapk_matrix  <- as.matrix(as_adjacency_matrix(g_mapk))
rownames(mapk_matrix) <- colnames(mapk_matrix)

prot_table2     <- cbind(prot_table,'site' = rep('',nrow(prot_table)))
prot_info_mapk  <- rbind(site_table[site_table$geneB %in% mapk_genes,],prot_table2[prot_table2$geneB %in% mapk_genes,])
prot_info_mapk2 <- unique(prot_info_mapk[prot_info_mapk$protA != prot_info_mapk$protB,])


write.csv(mapk_matrix,paste(mapk_result_path,out_adj_m,sep = '/'))
write.csv(prot_info_mapk2,paste(mapk_result_path,out_prot_table,sep = '/'),row.names = F)
write.csv(all_mapk2,paste(mapk_result_path,out_gene_table,sep = '/'),row.names = F)




