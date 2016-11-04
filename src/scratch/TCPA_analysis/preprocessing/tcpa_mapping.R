# Extract info from tcpa raw file
tcpa_mapping_raw    <- '~/Documents/workspace/phospho_network/RAWDATA/TCPA/TCPA_ANTIBODY_MAPPING.txt'
tcpa_mapping_outfile <- '~/Documents/workspace/phospho_network/processed_data/tcpa/tcpa_mapping.csv'

tcpa_mapping <- read.table(tcpa_mapping_raw,sep = '\t',as.is = T)

all_genes <- tcpa_mapping$V2
# correct errors in list:
all_genes[all_genes == 'Rab, 25'] <- 'RAB25'
all_genes2 <- gsub(',$','',all_genes)
all_genes3 <- gsub(',',';',all_genes)

all_positions <- lapply(strsplit(tcpa_mapping$V1,split = '_p'),function(x)x[-1])
all_positions2 <- lapply(all_positions,function(x)paste(x,collapse = '_'))
all_positions3 <- lapply(all_positions2,function(x)gsub(pattern = '[A-Za-z]',replacement = '',x))
all_positions4 <- lapply(all_positions3,function(x)gsub(pattern = '_',replacement = ';',x))
tcpa_info <- cbind('antibody' = tcpa_mapping$V1,'gene_symbol' = all_genes3,'site' = all_positions4)
write.csv(tcpa_info,tcpa_mapping_outfile,row.names = F,quote = T)
