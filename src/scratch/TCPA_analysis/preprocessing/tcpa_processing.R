# tcpa data processing
tcpa_info_file    <- '~/Documents/workspace/phospho_network/processed_data/tcpa/tcpa_mapping.csv'
tcpa_data_file    <- '~/scratch/TCPA_2016-03-22/TCGA-BRCA-L3-S35.csv'
tcpa_info_outfile <- '~/Documents/workspace/phospho_network/processed_data/tcpa/tcpa_data_processed.csv'

data_parser <- function(tcpa_data,tcpa_info){
  tcpa_data2 <- tcpa_data[-(1:2),]
  antibody <- tcpa_data[-(1:2),1]
  colnames(tcpa_data2) <- c('antibody',tcpa_data[1,-1])
  tcpa_data3 <- tcpa_data2[tcpa_data2[,1] %in% tcpa_info[,'antibody'],]
  return_matrix <- matrix(0,nrow = 0,ncol = 2+ncol(tcpa_data3))
  colnames(return_matrix) <- c(colnames(tcpa_info),colnames(tcpa_data3)[-1])
  for(i in 1:nrow(tcpa_data3)){
    data_info <- tcpa_info[tcpa_info[,1]==tcpa_data3[i,1],]
    genes <- strsplit(data_info$gene_symbol,split = ', ')[[1]]
    if(length(genes) == 1){
      return_matrix <- rbind(return_matrix,c(data_info,tcpa_data3[i,-1]))
    }else if(length(genes) > 1){
      data_info2 <- matrix(rep(data_info,length(genes)),nrow = length(genes),byrow = T)
      data_info2[,2] <- genes
      return_matrix <- rbind(return_matrix,cbind(data_info2,matrix(rep(tcpa_data3[i,-1],length(genes)),nrow = length(genes),byrow = T)))
    }
  }
  return(return_matrix)
}

#mannul correction:
data_correction <- function(tcpa_data_processed){
  tcpa_data_processed[tcpa_data_processed[,'gene_symbol'] == 'PIK3R1/2','gene_symbol'] <- 'PIK3R1;PIK3R2'
  tcpa_data_processed[,'gene_symbol'] <- gsub(';$','',tcpa_data_processed[,'gene_symbol'])
  tcpa_data_processed[,'gene_symbol'] <- gsub('^ ','',tcpa_data_processed[,'gene_symbol'])
  return(tcpa_data_processed)
}

#main body
tcpa_info           <- read.csv(tcpa_info_file,as.is = T)
tcpa_data_raw       <- t(read.csv(tcpa_data_file,header = F, as.is = T))
tcpa_data_processed <- data_parser(tcpa_data_raw,tcpa_info)
tcpa_data_corrected <- data_correction(tcpa_data_processed)
write.csv(tcpa_data_corrected,tcpa_info_outfile,row.names = F)
