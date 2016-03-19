# parse all KEGG pathways in human (hsa)
require(KEGGREST)
require(KEGGgraph)
require(reshape)
require(org.Hs.eg.db)
library(plyr)
kegg_raw_data_path = '~/Documents/workspace/phospho_network/RAWDATA/interaction_data/KEGG'
interaction_path   = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/processed_interaction_data'
kegg_out_file      = 'kegg_table_all.csv'

#helper function:
return_type <- function(x){
 if(is.null(x@subtype$subtype)){
   return('')
 }else{
   return(x@subtype$subtype@name)
 }
}

pathways  <- unique(unname(keggLink("pathway", "hsa")))
pathways2 <- gsub('path:(\\w+)','\\1',pathways)
interaction_table_kegg_raw <- NULL
for(i in pathways2){
  xml_file <- paste(paste(kegg_raw_data_path,i,sep = '/'),'xml',sep = '.')
  # download data, need this for first time run
  # retrieveKGML(i,organism = 'hsa',method = 'internal',destfile = xml_file)
  if(file.size(xml_file) > 1){
    pathway <- parseKGML2Graph(xml_file,expandGenes=TRUE)
    pathway_obj <- parseKGML(xml_file)
    pathway_name <- pathway_obj@pathwayInfo@title
    interactions <- getKEGGedgeData(pathway)
    interactions_table <- lapply(interactions,return_type)
    pathway_edges <- edges(pathway)
    edge_table <- melt(pathway_edges[lapply(pathway_edges,length)>0])
    if(!is.null(edge_table)){
      new_pathway_table <- cbind(edge_table,'pathway' = rep(pathway_name,nrow(edge_table)),'type' = rep('',nrow(edge_table)))
      new_pathway_table <- as.matrix(apply(new_pathway_table,2,as.character))
      if(ncol(new_pathway_table) == 1){
        new_pathway_table <- t(new_pathway_table)
      }
      for (j in 1:length(interactions_table)){
        interactors <- strsplit(names(interactions_table)[j],split = '~')[[1]]
          new_pathway_table[new_pathway_table[,'L1'] == interactors[1] & new_pathway_table[,'value'] == interactors[2],'type'] <- as.character(interactions_table[j])
        }
      interaction_table_kegg_raw <- rbind(interaction_table_kegg_raw,new_pathway_table)
    }
  }
}

#mapping kegg id to uniprot
interaction_table_kegg_raw2 <- apply(interaction_table_kegg_raw,2,as.character)

uniprot_mapping  <- keggConv("hsa","uniprot")
gene_mapping     <- keggConv("hsa","ncbi-geneid")

pkegg_ids       <- unname(uniprot_mapping)
gkegg_ids       <- unname(gene_mapping)

prot_ids        <- gsub('up:(\\w+)','\\1',names(uniprot_mapping))
gene_ids        <- gsub('ncbi-geneid:(\\w+)','\\1',names(gene_mapping))

uniprot_a <- mapvalues(interaction_table_kegg_raw2[,'L1'],from = pkegg_ids,to = prot_ids,warn_missing = F)
uniprot_b <- mapvalues(interaction_table_kegg_raw2[,'value'],from = pkegg_ids,to = prot_ids,warn_missing = F)
ncbi_a    <- mapvalues(interaction_table_kegg_raw2[,'L1'],from = gkegg_ids,to = gene_ids,warn_missing = F)
ncbi_b    <- mapvalues(interaction_table_kegg_raw2[,'value'],from = gkegg_ids,to = gene_ids,warn_missing = F)
uniprot_a2<- gsub('hsa:\\w+','',uniprot_a)
uniprot_b2<- gsub('hsa:\\w+','',uniprot_b)

gene_mapping_a <- select(org.Hs.eg.db, keys=ncbi_a, columns=c('SYMBOL'), keytype="ENTREZID")
gene_mapping_b <- select(org.Hs.eg.db, keys=ncbi_b, columns=c('SYMBOL'), keytype="ENTREZID")

interaction_table_kegg <- cbind(gene_mapping_a$SYMBOL,uniprot_a2,gene_mapping_b$SYMBOL,uniprot_b2,interaction_table_kegg_raw2[,'pathway'],interaction_table_kegg_raw2[,'type'])
colnames(interaction_table_kegg) <- c('geneA','protA','geneB','protB','pathway','type')
interaction_table_kegg <- unique(apply(interaction_table_kegg,2,as.character))
write.csv(interaction_table_kegg,paste(interaction_path,kegg_out_file,sep = '/'),row.names = F)
