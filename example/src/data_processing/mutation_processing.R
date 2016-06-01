# missense mutation occured more than missense_occurance or occured in mutation_effect_file would be treated as an independent predictor
# other types of mutations would be treated as a binary variable "mutation_misc"
mutation_file         = '~/Documents/workspace/phospho_network/example/RAWDATA/brca_mutation_raw.maf' #maf format
mutation_effect_file  = '~/Documents/workspace/phospho_network/example/RAWDATA/mutations_effect.tsv' #tab delimited file with "Gene	Substitution	Effect"
mutation_misc_out     = '~/Documents/workspace/phospho_network/example/script_files/mutation_matrix_misc.csv'
mutation_missense_out = '~/Documents/workspace/phospho_network/example/script_files/mutation_missense.csv'
rppa_file             = '~/Documents/workspace/phospho_network/example/script_files/rppa_processed.csv'
missense_occurance    = 1 

mut_data_raw   <- read.table(mutation_file,sep = '\t',skip = 1,as.is = T,header = F,fill = NA)
mut_effect_raw <- read.table(mutation_effect_file,sep = '\t',header = T)
sample_rppa    <- colnames(read.csv(rppa_file,row.names = 1))

mut_data <- mut_data_raw[-1,1:ncol(mut_data_raw)-1]
colnames(mut_data) <- mut_data_raw[1,1:ncol(mut_data_raw)-1]
mut_data <- mut_data[grep(pattern = '^LOC\\d+',x = ,mut_data$Hugo_Symbol,invert = T),]
mut_data <- mut_data[grep(pattern = '^ENSG\\d+',x = ,mut_data$Hugo_Symbol,invert = T),]
mut_data$Tumor_Sample_Barcode <- gsub(pattern = '(\\w+-\\w+-\\w+-\\d+).+',replacement = '\\1',x = mut_data$Tumor_Sample_Barcode)
mut_data <- mut_data[mut_data$Tumor_Sample_Barcode %in% gsub('\\.','-',sample_rppa),]

mutation_missense <- mut_data[mut_data$Variant_Classification == 'Missense_Mutation',]
mutation_misc     <- mut_data[!mut_data$Variant_Classification == 'Missense_Mutation',]
mut_effect <- mut_effect_raw[mut_effect_raw$Effect != 'nt',]

gene_list <- unique(mutation_misc$Hugo_Symbol)
sample_list <- unique(mut_data$Tumor_Sample_Barcode)
mutation_matrix_misc <- matrix(0,nrow = length(gene_list),ncol = length(sample_list))
colnames(mutation_matrix_misc) <- sample_list
rownames(mutation_matrix_misc) <- gene_list
for (i in 1:nrow(mutation_misc)){
  mutation_matrix_misc[mutation_misc$Hugo_Symbol[i],mutation_misc$Tumor_Sample_Barcode[i]] <- 1
}

missense_names  <- paste(mutation_missense$Hugo_Symbol,gsub('p\\.','\\1',mutation_missense$amino_acid_change_WU),sep = '_')
missense_table  <- table(missense_names)
missense_keep   <- unique(names(missense_table[missense_table > missense_occurance]),intersect(names(missense_table),paste(mut_effect$Gene,mut_effect$Substitution,sep = '_')))
missense_matrix <- cbind('gene' = gsub('(\\w+)_\\w+','\\1',missense_names),missense_names,'sample'=mutation_missense$Tumor_Sample_Barcode)[missense_names %in% missense_keep,]

write.csv(mutation_matrix_misc,mutation_misc_out)
write.csv(missense_matrix,mutation_missense_out)


