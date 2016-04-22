# transform TCGA level 3 CNV data to gene and mean 
cnv_file = '~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/brca_hg19_qc.merged.seg'
rna_file = '~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/BRCA.exp.547.med.txt'

library("biomaRt")
cnv_data <- read.table(cnv_file,as.is = T,header = T,sep = '\t')
rna_data <- read.delim(rna_file,row.names = 1)


rna_samples <- gsub('\\.','-',colnames(rna_data))
rna_samples2 <- gsub('(\\w+-\\w+-\\w+-\\d+).+','\\1',rna_samples)
cnv_samples <- unique(cnv_data$Sample)
colnames(rna_data) <- rna_samples2

intersect_samples <- intersect(rna_samples2,cnv_samples)
cnv_data_intersect_raw <- cnv_data[cnv_data$Sample %in% intersect_samples,]
rna_data_intersect <- rna_data[,intersect_samples]

ensembl54=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

return_gene <- function(x){
  results=getBM(attributes = c("hgnc_symbol"),
                filters = c("chromosomal_region","biotype"),values = list(paste(x[2:4],collapse = ':'),"protein_coding"), mart = ensembl54)
  genes <- paste(results$hgnc_symbol,collapse = ';')
  return(genes)
}

cnv_data_intersect <- cnv_data_intersect_raw[1:5,]

filterlist <- apply(cnv_data_intersect,1,return_gene)
cnv_mapping <- cbind.data.frame(cnv_data_intersect,filterlist, stringsAsFactors = FALSE)

cnv_gene <- matrix(0,ncol = 3,nrow = 0)
samples <- unique(cnv_mapping$Sample)
for (i in samples){
  m = cnv_mapping[cnv_mapping$Sample == i,]
  new_table <- matrix(0,nrow = 0,ncol = 3)
  for (j in 1:nrow(m)){
    new_genes <- strsplit(m[1,]$filterlist,split = ';')[[1]]
    if(length(new_genes) != 0){
      new_table <- rbind(new_table,cbind.data.frame('sample' = i,'gene' = new_genes,'seg.mean' = m[1,]$seg.mean))
    }
  }
  if(nrow(new_table) > 0){
    cnv_gene <- rbind(cnv_gene,aggregate(seg.mean ~ sample + gene, FUN = "mean", data = new_table)) #mean
  }
}
cnv_gene_matrix <- matrix(NA,nrow = length(unique(cnv_gene$sample)),ncol = length(unique(cnv_gene$gene)))
colnames(cnv_gene_matrix) <- unique(cnv_gene$gene)
rownames(cnv_gene_matrix) <- unique(cnv_gene$sample)

for (i in 1: nrow(cnv_gene)){
  new_line <- cnv_gene[i,]
  cnv_gene_matrix[new_line$sample,new_line$gene] <- as.numeric(new_line$seg.mean)
}
cnv_gene_matrix <- t(cnv_gene_matrix)
genes_intersect <- intersect(rownames(cnv_gene_matrix),rownames(rna_data_intersect))
sample_intersect <- intersect(colnames(cnv_gene_matrix),colnames(rna_data_intersect))
rna_matrix_intersect <- rna_data_intersect[genes_intersect,sample_intersect]
cnv_matrix_intersect <- cnv_gene_matrix[genes_intersect,sample_intersect]

write.csv(rna_matrix_intersect,'rna_processed.csv')
write.csv(cnv_matrix_intersect,'cnv_processed.csv')


#with all genes
# cors <- c()
# max_abs <- c()
# iqrs <- c()
# for (i in 1:nrow(rna_matrix_intersect)){
#   ind <- (!is.na(rna_matrix_intersect[i,])) & (!is.na(cnv_matrix_intersect[i,]))
#   cors <- c(cors,cor(as.numeric(rna_matrix_intersect[i,ind]),as.numeric(cnv_matrix_intersect[i,ind]),method = 'spearman'))
#   max_abs <- c(max_abs,max(abs(as.numeric(cnv_matrix_intersect[i,ind]))))
#   iqrs <- c(iqrs,IQR(as.numeric(cnv_matrix_intersect[i,ind])))
# }

