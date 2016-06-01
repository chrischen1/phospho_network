measure_filename      = '~/Documents/workspace/phospho_network/RAWDATA/standard_example/brca_RPPA_raw.txt'
rppa_mapping_file1    = '~/Documents/workspace/phospho_network/RAWDATA/standard_example/brca_rppa_mapping1.txt'
rppa_mapping_file2    = '~/Documents/workspace/phospho_network/RAWDATA/standard_example/brca_rppa_mapping2.txt'

outfile_name          = '~/Documents/workspace/phospho_network/script_files/standard_example/rppa_processed_raw.csv'

rppa_mapping <- function(measure_filename,rppa_mapping_file1,rppa_mapping_file2){
  mdata <- read.table(measure_filename,sep = '\t',as.is = T)
  mappping1 <- read.table(rppa_mapping_file1,sep = '\t',as.is = T)
  mappping2 <- read.table(rppa_mapping_file2,sep = '\t',as.is = T)
  mdata2 <- mdata[-1,-1]
  colnames(mdata2) <- gsub(pattern = '(\\w+\\.\\w+\\.\\w+\\.\\d+).+',replacement = '\\1',x = mdata[1,-1])
  #mapping gene names
  antibody_names <- tolower(gsub('(.+)-[A-Z]+-[A-Z]+','\\1',mappping1$V3[-1]))
  antibody_table <- rbind(cbind(mappping1$V1[-1],antibody_names),cbind(mappping2$V2,tolower(mappping2$V1)))
  antibody_table[,2] <- gsub('[\\.|\\-]','_',antibody_table[,2])
  antibody_table <- unique(antibody_table)
  antibody_list  <- gsub('[\\.|\\-]','_',tolower(mdata[-1,1]))
  #mapping position
  positions <- lapply(strsplit(antibody_table[,2],split = '_ps|_pt|_py'),function(x)gsub('[a-z]','',paste(x[-1],collapse = '_')))
  positions2 <- gsub('_$','',unlist(positions),perl = T)
  # correct errors in list:
  antibody_table[antibody_table == 'Rab, 25'] <- 'RAB25'
  genes <- gsub('\\W+',';',antibody_table[,1])
  genes[genes == 'PIK3R1;2'] <- 'PIK3R1;PIK3R2'
  genes2 <- gsub(';$','',genes)
  #paste gene name and position
  gene_pos <- paste(genes2,gsub('_',';',positions2),sep = '_')
  names(gene_pos) <- antibody_table[,2]
  gene_pos['bcl_xl'] <- "obsolutedBCL2L1_"
  rownames(mdata2) <- gene_pos[antibody_list]
  return(mdata2)
}

mdata <- rppa_mapping(measure_filename,rppa_mapping_file1,rppa_mapping_file2)
write.csv(mdata,outfile_name)
