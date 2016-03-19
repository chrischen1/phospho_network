# analysis the coverage of interactions from different data sources:
library(gplots)
interaction_path  = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed'
processed_path    = paste(interaction_path,'processed_interaction_data',sep = '/')

ms_gene_list = 'gene_experiment_list.txt'
signor_file  = 'signor_table_all.csv'
ks_file      = 'ks_table_all.csv'
bg_file      = 'bg_table_all.csv'
pw_file      = 'pw_table_all.csv'

plot_venn <- function(signor_data,ks_data,bg_data,pw_data){
  col_nameA   <- colnames(bg_data)[1]
  col_nameB   <- colnames(bg_data)[2]
  int_signor  <- paste(signor_data[,col_nameA],signor_data[,col_nameB],sep = '_')
  int_ks      <- paste(ks_data[,col_nameA],ks_data[,col_nameB],sep = '_')
  int_bg      <- paste(bg_data[,1],bg_data[,2],sep = '_')
  int_pw      <- paste(pw_data[,1],bg_data[,2],sep = '_')
  v.table<-venn(list(signor=int_signor,kinase_substrate=int_ks,biogrid = int_bg,pathway_commons=int_pw),show.plot = F)
}

par(mfrow=c(2,1))

ms_genes    <- unlist(read.table(paste(interaction_path,ms_gene_list,sep = '/'),as.is = T))
signor_data <- read.csv(paste(processed_path,signor_file,sep = '/'),as.is = T)
ks_data     <- read.csv(paste(processed_path,ks_file,sep = '/'),as.is = T)
bg_data     <- read.csv(paste(processed_path,bg_file,sep = '/'),as.is = T)
pw_data     <- read.csv(paste(processed_path,pw_file,sep = '/'),as.is = T)

signor_data_ms <- signor_data[signor_data$geneA %in% ms_genes & signor_data$geneB %in% ms_genes,]
ks_data_ms     <- ks_data[ks_data$geneA %in% ms_genes & ks_data$geneB %in% ms_genes,]
bg_data_ms     <- bg_data[bg_data$geneA %in% ms_genes & bg_data$geneB %in% ms_genes,]
pw_data_ms     <- pw_data[pw_data$geneA %in% ms_genes & pw_data$geneB %in% ms_genes,]

ms_plot  <- plot_venn(signor_data_ms,ks_data_ms,bg_data_ms,pw_data_ms)
all_plot <- plot_venn(signor_data,ks_data,bg_data,pw_data)

plot(ms_plot)
text(210,400,'coverage of interactions in MS data')
plot(all_plot)
text(210,400,'coverage of interactions of all data')
