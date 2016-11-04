library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#read info about go term that need retrieving
act <- read.table('~/Documents/workspace/phospho_network/RAWDATA/TF_go/tf_activator_terms.txt',sep = '\t',as.is = T)
rep <- read.table('~/Documents/workspace/phospho_network/RAWDATA/TF_go/tf_repressor_terms.txt',sep = '\t',as.is = T)


results <- getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"), filters = "go_id",values = act$V1, mart = mart,uniqueRows = T)
results2 <- getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"), filters = "go_id",values = rep$V1, mart = mart,uniqueRows = T)


ambi <- intersect(results$hgnc_symbol,results2$hgnc_symbol)
act_table <- cbind(results[!results$hgnc_symbol %in% ambi,], 'effect' = 'activator')
rep_table <- cbind(results[!results2$hgnc_symbol %in% ambi,], 'effect' = 'repressor')
final_table <- rbind(act_table,rep_table)

write.csv(final_table,'~/Documents/workspace/phospho_network/script_files/TF_analysis/GO_tf_table.csv')
