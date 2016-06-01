# parse all KEGG pathways in human (hsa)
require(KEGGREST)
require(KEGGgraph)
require(reshape)
require(org.Hs.eg.db)
library(plyr)
kegg_raw_data_path = '~/Documents/workspace/phospho_network/example/script_files/network/KEGG'
kegg_out_path      = '~/Documents/workspace/phospho_network/example/script_files/network'
kegg_out_file      = 'kegg_table_all.csv'

pathways  <- unique(unname(keggLink("pathway", "hsa")))
pathways2 <- gsub('path:(\\w+)','\\1',pathways)
interaction_table_kegg_raw <- NULL
for(i in pathways2){
  xml_file <- paste(paste(kegg_raw_data_path,i,sep = '/'),'xml',sep = '.')
  retrieveKGML(i,organism = 'hsa',method = 'internal',destfile = xml_file) # download data, need this for first time run
  if(file.size(xml_file) > 1){
    pathway <- parseKGML2DataFrame(xml_file)
    pathway_obj <- parseKGML(xml_file)
    pathway_name <- pathway_obj@pathwayInfo@title
    interaction_table_kegg_raw <- rbind(interaction_table_kegg_raw,cbind(pathway,'pathway' = rep(pathway_name,nrow(pathway))))
  }
}
interaction_table_kegg_raw2 <- data.frame(apply(interaction_table_kegg_raw,2,as.character),stringsAsFactors = F)

#mapping kegg id to uniprot
uniprot_mapping  <- keggConv("hsa","uniprot")
gene_mapping     <- keggConv("hsa","ncbi-geneid")

pkegg_ids       <- unname(uniprot_mapping)
gkegg_ids       <- unname(gene_mapping)

prot_ids        <- gsub('up:(\\w+)','\\1',names(uniprot_mapping))
gene_ids        <- gsub('ncbi-geneid:(\\w+)','\\1',names(gene_mapping))

uniprot_a <- mapvalues(interaction_table_kegg_raw2$from,from = pkegg_ids,to = prot_ids,warn_missing = F)
uniprot_b <- mapvalues(interaction_table_kegg_raw2$to,from = pkegg_ids,to = prot_ids,warn_missing = F)
ncbi_a    <- mapvalues(interaction_table_kegg_raw2$from,from = gkegg_ids,to = gene_ids,warn_missing = F)
ncbi_b    <- mapvalues(interaction_table_kegg_raw2$to,from = gkegg_ids,to = gene_ids,warn_missing = F)
uniprot_a2<- gsub('hsa:\\w+','',uniprot_a)
uniprot_b2<- gsub('hsa:\\w+','',uniprot_b)

gene_mapping_a <- select(org.Hs.eg.db, keys=ncbi_a, columns=c('SYMBOL'), keytype="ENTREZID")
gene_mapping_b <- select(org.Hs.eg.db, keys=ncbi_b, columns=c('SYMBOL'), keytype="ENTREZID")


interaction_table_kegg <- cbind(gene_mapping_a$SYMBOL,uniprot_a2,gene_mapping_b$SYMBOL,uniprot_b2,interaction_table_kegg_raw2$pathway,interaction_table_kegg_raw2$subtype)
colnames(interaction_table_kegg) <- c('geneA','protA','geneB','protB','pathway','type')
interaction_table_kegg <- unique(apply(interaction_table_kegg,2,as.character))
write.csv(interaction_table_kegg,paste(kegg_out_path,kegg_out_file,sep = '/'),row.names = F)
