# msdata_processing by RPPA format

msdata_file <- '~/Documents/workspace/phospho_network/script_files/msdata_processed.csv'
mstot_file  <- '~/Documents/workspace/phospho_network/script_files/total_protein_processed.csv'
rppa_file   <- '~/Documents/workspace/phospho_network/script_files/analysis_rppa/result_matrix.csv'

overlap_out <- '~/Documents/workspace/phospho_network/script_files/analysis_rppa/response_list.txt'

msdata  <- read.csv(msdata_file,as.is = T)
totdata <- read.csv(mstot_file,as.is = T)

rppa_data    <- read.csv(rppa_file,as.is = T)
msdata_yname <- paste(msdata$gene_symbol,msdata$Site.Position,sep = '_')
rppa_yname   <- unique(rppa_data$gene_site)
rppa_gene    <- gsub('(\\w+)_.+','\\1',rppa_yname)

overlap_y    <- intersect(msdata_yname,rppa_yname)




write.table(overlap_y,overlap_out,row.names = F,col.names = F)
