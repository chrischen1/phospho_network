vepDir  = '/groups/sorger/cchris/vep_files'
outfile = '/groups/sorger/cchris/cancer_gene_mutation.tsv'
gene_list_file = '~/Documents/workspace/phospho_network/RNAseq/cellines_info/CancerGenes.tsv'
sample_info_file = '~/Documents/workspace/phospho_network/RNAseq/annotation/RNAseq_sample.tsv'
variant_keep = c('coding_sequence_variant,NMD_transcript_variant','frameshift_variant','frameshift_variant,NMD_transcript_variant',
                 'incomplete_terminal_codon_variant,coding_sequence_variant','inframe_deletion','inframe_deletion,NMD_transcript_variant',
                 'inframe_insertion','inframe_insertion,NMD_transcript_variant','missense_variant','missense_variant,NMD_transcript_variant',
                 'missense_variant,splice_region_variant','missense_variant,splice_region_variant,NMD_transcript_variant','start_lost ')

gene_list <- read.table(gene_list_file,header = T,sep = '\t',as.is = T)
sample_info <- read.table(sample_info_file,as.is = T,sep = '\t',header = T)
sample_info2 <- sample_info$cellines
names(sample_info2) <- sample_info$sample_id
  
gene_names <- gene_list$Gene
vep_files <- grep('.vep.txt$',list.files(vepDir),value = T)
all_table <- NULL
for(i in vep_files){
  sample_id = gsub('.vep.txt','',i,fixed = T)
  vep_table <- read.table(paste(vepDir,i,sep = '/'),as.is = T)
  gene_symbol <- gsub('.+SYMBOL\\=(.+);SYMBOL_SOURCE.+','\\1',vep_table[,14])
  vep_table <- cbind(gene_symbol,vep_table,'sample_id' = sample_id)
  gene_names_valid <- intersect(gene_names,gene_symbol)
  for (gene in gene_names_valid){
    all_table <- rbind(all_table,vep_table[vep_table$gene_symbol==gene,])
  }
}
result_table <- all_table
result_table$sample_id <- sample_info2[gsub(pattern = '_QC','',toupper(all_table$sample_id))]
result_table$sample_id[is.na(result_table$sample_id)] <- all_table$sample_id[is.na(result_table$sample_id)]
cancer_gene_mutation <- result_table[result_table$V7 %in% variant_keep,c(1,8,9,10,11,12,ncol(result_table))]
write.table(cancer_gene_mutation,outfile,sep = '\t')
