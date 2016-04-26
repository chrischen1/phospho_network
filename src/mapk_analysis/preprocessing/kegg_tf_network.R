library(biomaRt)
kegg_table <- read.csv('~/Documents/workspace/phospho_network/script_files/interaction_data_processed/processed_interaction_data/kegg_table_all.csv')
TF_KEGG    <- read.table('~/Documents/workspace/phospho_network/RAWDATA/TF_kegg/hsa03000.keg',sep = '\n')
outfile    <- '~/Documents/workspace/phospho_network/script_files/TF_analysis/tf_table.csv'
types      <- c('activation','inhibition','phosphorylation','dephosphorylation','repression')
pactive    <- c('activation_phosphorylation','dephosphorylation_inhibition')
pinhibit   <- c('activation_dephosphorylation','inhibition_phosphorylation')

tf_table <- matrix(0,nrow = 0, ncol = 3)
colnames(tf_table) <- c('TF','active_state','phenotype')

TF_list_raw <- gsub(pattern = '^D\\s+\\d+\\s+(.+); ',replacement = '\\1',grep(pattern = '^D\\s+',TF_KEGG$V1,value = T))
TF_list     <- unique(unlist(lapply(strsplit(TF_list_raw,split = ';'),function(x)x[1])))

kegg_table_tf  <- unique(kegg_table[kegg_table$geneB %in% TF_list & kegg_table$type %in% types, c('geneA','geneB','type')])
kegg_table_tf$type <- gsub('repression','inhibition',kegg_table_tf$type)

kegg_table_tf2 <- aggregate(type ~ geneA + geneB, FUN = function(x)paste(sort(x),collapse = '_'), data = kegg_table_tf)
kegg_table_tf3 <- kegg_table_tf2[kegg_table_tf2$type %in% c(pactive,pinhibit),]

TF_list_filtered <- unique(as.character(kegg_table_tf3$geneB))

for (tf in TF_list_filtered){
  npa <- sum((kegg_table_tf3$type[kegg_table_tf3$geneB == tf]) %in% pactive)
  npi <- sum(kegg_table_tf3$geneB == tf) - npa
  # vote
  if(npa > npi){
    tf_table <- rbind(tf_table,c(tf,'p',''))
  }else if(npa < npi){
    tf_table <- rbind(tf_table,c(tf,'np',''))
  }
}
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
for (i in 1:nrow(tf_table)){
  results <- getBM(attributes = c("go_id","name_1006",'namespace_1003'), filters = "hgnc_symbol",values = tf_table[i,'TF'], mart = mart)
  
  tf_table[i,'phenotype'] <- paste(results$name_1006[results$namespace_1003 == 'molecular_function'],collapse = ';')
}

write.csv(tf_table,  outfile)
  
  
  
  