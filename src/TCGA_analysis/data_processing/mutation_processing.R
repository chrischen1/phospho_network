mutation_file = '~/Documents/workspace/phospho_network/RAWDATA/tcga_brac/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf'

mut_data_raw <- read.table(mutation_file,sep = '\t',skip = 1,as.is = T,header = F,fill = NA)
info_ind <- c(1,2,5,8,9,10,16,17)
mut_data <- mut_data_raw[-1,info_ind]
colnames(mut_data) <- mut_data_raw[1,info_ind]

gene_list <- unique(mut_data$Hugo_Symbol)
sample_list <- unique(mut_data$Tumor_Sample_Barcode)

mutation_matrix <- matrix('',nrow = length(gene_list),ncol = length(sample_list))
colnames(mutation_matrix) <- sample_list
rownames(mutation_matrix) <- gene_list

for (i in 1:nrow(mut_data)){
  mutation_matrix[mut_data$Hugo_Symbol[i],mut_data$Tumor_Sample_Barcode[i]] <- paste(mutation_matrix[mut_data$Hugo_Symbol[i],mut_data$Tumor_Sample_Barcode[i]],mut_data$Variant_Classification[i],';',sep = '')
}

write.csv(mutation_matrix,'~/Documents/workspace/phospho_network/script_files/TCGA_analysis/mutation_matrix_classified.csv')


