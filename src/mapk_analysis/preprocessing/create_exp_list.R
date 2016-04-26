script_file_path     = '/Users/xxx/Documents/workspace/phospho_network/script_files'
interactions_path    = '/Users/xxx/Documents/workspace/phospho_network/script_files/interaction_data_processed'

total_protein_file   = 'total_protein_processed.csv'
ms_protein_file      = 'msdata_processed.csv'

tot_data <- read.csv(paste(script_file_path,total_protein_file,sep='/'),as.is = T)
ms_data  <- read.csv(paste(script_file_path,ms_protein_file,sep='/'),as.is = T)

intersect_protein <- intersect(unique(tot_data$merge_id),unique(ms_data$merge_id))
intersect_gene    <- intersect(unique(tot_data$Gene.Symbol),unique(ms_data$gene_symbol))
intersect_gene2   <- intersect_gene[intersect_gene != '']
write.table(intersect_protein,file = paste(interactions_path,'prot_experiment_list.txt',sep='/'),quote = F,row.names = F,col.names = F)
write.table(intersect_gene2,file = paste(interactions_path,'gene_experiment_list.txt',sep='/'),quote = F,row.names = F,col.names = F)
