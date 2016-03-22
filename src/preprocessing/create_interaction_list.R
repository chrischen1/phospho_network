# This script would analysis interaction information from different database for all proteins/genes, 
# and return 2 long table which rows are [kinase substrate site] or [kinase substrate].
# Autointeractions are included!!!
interaction_data_path = '~/Documents/workspace/phospho_network/RAWDATA/interaction_data'
interactions_path     = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed'
processed_path        = paste(interactions_path,'processed_interaction_data',sep = '/')

ks_data_file          = 'Kinase_Substrate_Dataset.csv'
signor_file           = 'all_data_SIGNOR.tsv'
bg_data_file          = 'BIOGRID-ALL-3.4.134.tab2.txt'
pw_data_file          = 'Pathway Commons.7.All.BINARY_SIF.hgnc.sif'
output_file_site      = 'interaction_table_site_all.csv'
output_file_prot      = 'interaction_table_prot_all.csv'
output_file_gene      = 'interaction_table_gene_all.csv'
output_file_network   = 'network_all_gene.csv'
kegg_table_file       = 'kegg_table_all.csv'

only_phos             = F #Include phosphorylation only?
PC_interactions       = c('interacts-with','controls-state-change-of','controls-phosphorylation-of')
kegg_interactions     = c('activation','dephosphorylation','inhibition','phosphorylation','repression','state change','ubiquitination')

# parse from signor:
signor_parser <-function(signor){
  signor_data <- read.table(signor,as.is = T,header = T,sep = '\t',fill = T,quote = "")
  signor_data_filtered <- signor_data[signor_data$TypeA=='PROTEIN' & signor_data$TypeB=='PROTEIN',c('EntityA','IdA','EntityB','IdB','Effect','EffectMechanism','MechanismResidues')]
  if(only_phos){
    signor_data_filtered <- signor_data_filtered[signor_data_filtered$EffectMechanism=='phosphorylation',]
  }
  signor_table <- cbind(signor_data_filtered$EntityA,signor_data_filtered$IdA,signor_data_filtered$EntityB,signor_data_filtered$IdB,gsub(pattern = '[A-Za-z]+(\\d+)',replacement = '\\1',signor_data_filtered$MechanismResidues))
  signor_table_unique <- unique(signor_table)
  return(signor_table_unique)
}


#parse from Kinase_Substrate_Dataset
ks_parser <- function(ks){
  ks_data <- read.csv(ks,header = T,as.is = T,skip = 2)
  ks_data_filtered <- ks_data[ks_data$SUB_ORGANISM == 'human' & ks_data$KIN_ORGANISM == 'human',]
  ks_data_table <- ks_data_filtered[,c('GENE','KIN_ACC_ID','SUB_GENE','SUB_ACC_ID','SUB_MOD_RSD')]
  ks_data_table$SUB_MOD_RSD <- gsub(pattern = '[A-Za-z]+',replacement = '',x = ks_data_table$SUB_MOD_RSD)
  ks_data_table_unique <- unique(ks_data_table)
  return(ks_data_table_unique)
}


#parse from biogrid
bg_parser <- function(bg){
  bg_data <- read.table(bg,sep = '\t',fill = T,as.is = T,header = T,quote = '')
  bg_data_filtered <- bg_data[bg_data$Organism.Interactor.A=='9606' & bg_data$Organism.Interactor.B=='9606',c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
  bg_data_table_unique <- unique(bg_data_filtered)
  return(bg_data_table_unique)
}

# parse from pathway commons
pw_parser <- function(pw){
  pw_data <- read.table(pw,fill = T,as.is = T,header = F,sep = '\t')
  pw_data_table <- pw_data[pw_data$V2 %in% PC_interactions,]
  pw_data_table_unique <- unique(pw_data_table[,c('V1','V3')])
  return(pw_data_table_unique)
}

#helper functions
# Cartesian_product for 2 genes by their uniprot accessions
cartesian_product <- function(x,y,ref_list){
  x_uniprot <- ref_list$To[ref_list$From==x]
  y_uniprot <- ref_list$To[ref_list$From==y]
  xy_table <- expand.grid(x_uniprot,y_uniprot)
  return(xy_table)
}

expand_matrix <- function(m,ncol,sep = ';'){
  m2 <- m[0,]
  for(i in 1:nrow(m)){
    var_split <- strsplit(m[i,ncol],split = ';')[[1]]
    if(length(var_split) <= 1){
      m2 <- rbind(m2,m[i,])
    }else if(length(var_split) > 1){
      new_table <- matrix(rep(as.character(m[i,]),length(var_split)),nrow = length(var_split),byrow = T)
      new_table[,ncol] <- var_split
      colnames(new_table) <- colnames(m2)
      m2 <- rbind(m2,new_table)
    }
  }
  return(m2)
}

# Main body
# Create table for interactions
interaction_table_site <- matrix(0,ncol=5,nrow = 0)
colnames(interaction_table_site) <- c('geneA','protA','geneB','protB','site')
interaction_table_prot <- matrix(0,ncol=4,nrow = 0)
colnames(interaction_table_prot) <- c('geneA','protA','geneB','protB')
interaction_table_gene <- matrix(0,ncol=2,nrow = 0)
colnames(interaction_table_gene) <- c('geneA','geneB')

# SIGNOR
if(signor_file != ''){
  signor_table <- signor_parser(paste(interaction_data_path,signor_file,sep = '/'))
  colnames(signor_table) <- colnames(interaction_table_site)  
  write.csv(signor_table,paste(processed_path,'signor_table_all.csv',sep = '/'),row.names = F)
}else{
  signor_table <- NULL
}

#Kinase-substrate
if(ks_data_file != ''){
  ks_table <- ks_parser(paste(interaction_data_path,ks_data_file,sep = '/'))
  colnames(ks_table) <- colnames(interaction_table_site)  
  write.csv(ks_table,paste(processed_path,'ks_table_all.csv',sep = '/'),row.names = F)
}else{
  ks_table <- NULL
}

#BioGrid
if(bg_data_file != ''){
  bg_table <- bg_parser(paste(interaction_data_path,bg_data_file,sep = '/'))
  colnames(bg_table) <- colnames(interaction_table_gene)
  write.csv(bg_table,paste(processed_path,'bg_table_all.csv',sep = '/'),row.names = F)
}else{
  bg_table <- NULL
}

#Pathway Commons
if(pw_data_file != ''){
  pw_table <- pw_parser(paste(interaction_data_path,pw_data_file,sep = '/'))
  colnames(pw_table) <- colnames(interaction_table_gene)
  write.csv(pw_table,paste(processed_path,'pw_table_all.csv',sep = '/'),row.names = F)
}else{
  pw_table <- NULL
}

# Combine data from different sources
interaction_table_site <- unique(rbind(signor_table[signor_table[,'site'] != '',],ks_table))
interaction_table_prot <- signor_table[signor_table[,'site'] == '',1:4]
interaction_table_gene <- unique(rbind(signor_table[signor_table[,'site'] == '',c(1,3)],bg_table,pw_table))
output_network         <- unique(rbind(interaction_table_site[,c('geneA','geneB')],interaction_table_prot[,c('geneA','geneB')],interaction_table_gene))
output_network2        <- output_network[output_network$geneA != '' & output_network$geneB != '',]

# Add info from KEGG
if(kegg_table_file != ''){
  kegg_table  <- read.csv(paste(processed_path,kegg_table_file,sep = '/'),as.is = T)
  kegg_table2 <- kegg_table[kegg_table$type %in% kegg_interactions, ]
  interaction_table_prot <- unique(rbind(interaction_table_prot,kegg_table2[,colnames(interaction_table_prot)]))
  interaction_table_gene <- unique(rbind(interaction_table_gene,kegg_table2[,colnames(interaction_table_gene)]))
  output_network2 <- unique(rbind(output_network2,kegg_table2[,colnames(output_network2)]))
}

interaction_table_site2 <- unique(expand_matrix(interaction_table_site,ncol = which(colnames(interaction_table_site)=='site')))

write.csv(interaction_table_site2,paste(interactions_path,output_file_site,sep = '/'),row.names = F)
write.csv(interaction_table_prot,paste(interactions_path,output_file_prot,sep = '/'),row.names = F)
write.csv(interaction_table_gene,paste(interactions_path,output_file_gene,sep = '/'),row.names = F)
write.csv(output_network2,paste(interactions_path,output_file_network,sep = '/'),row.names = F)





