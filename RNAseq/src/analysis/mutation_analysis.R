RNAseq_sample <- read.delim("~/Documents/workspace/phospho_network/RNAseq/annotation/RNAseq_sample.tsv", stringsAsFactors=FALSE)
result_table <- read.delim("~/Documents/workspace/phospho_network/RNAseq/processed_data/result_table.txt", stringsAsFactors=FALSE)
PIK3CA_mutations <- read.delim("~/Documents/workspace/phospho_network/RNAseq/annotation/PIK3CA_mutations.txt")

samples <- RNAseq_sample$cellines
names(samples) <- RNAseq_sample$sample_id
result_table$sample_id <- samples[gsub('_qc','',result_table$sample_id)]
mutation_table <- result_table[result_table$V7=='missense_variant',]

mutation_ref <- PIK3CA_mutations[grep('BREAST',PIK3CA_mutations$Sample),]
mutation_ref$Sample <- gsub(' BREAST','',mutation_ref$Sample)
mutation_ref_exp <- mutation_ref[mutation_ref$Sample %in% RNAseq_sample$cellines,]

write.table(mutation_ref_exp,'~/Documents/workspace/phospho_network/RNAseq/processed_data/mutation_ref.tsv',sep = '\t',quote = F,row.names = F)
write.table(mutation_table,'~/Documents/workspace/phospho_network/RNAseq/processed_data/mutation_exp.tsv',sep = '\t',quote = F,row.names = F,col.names = F)
