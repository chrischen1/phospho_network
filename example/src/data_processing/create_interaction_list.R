# This script would analysis interaction information from different database for all proteins/genes, 
# and return 2 long table which rows are [kinase substrate site] or [kinase substrate].
interaction_data_path = '~/Documents/workspace/phospho_network/example/RAWDATA/interactions'
output_path           = '~/Documents/workspace/phospho_network/example/script_files/network'

ks_data_file          = 'Kinase_Substrate_Dataset.csv'
signor_file           = 'all_data_SIGNOR.tsv'
bg_data_file          = ''      #BIOGRID-ALL-3.4.134.tab2.txt if PPIs are considered
pw_data_file          = 'Pathway Commons.7.All.BINARY_SIF.hgnc.sif'
output_file_site      = 'interaction_table_site_all.csv'
output_file_prot      = 'interaction_table_prot_all.csv'
output_file_gene      = 'interaction_table_gene_all.csv'
output_file_network   = 'network_all_gene.csv'
kegg_table_file       = 'kegg_table_all.csv'

signor_interactions   = c('dephosphorylation','phosphorylation')
PC_interactions       = c('in-complex-with','controls-state-change-of','controls-phosphorylation-of')
kegg_interactions     = c('activation','dephosphorylation','inhibition','phosphorylation','repression','state change','compound')
include_auto          = F #include autointeractions?

# parse from signor:
signor_parser <-function(signor,signor_interactions,include_auto){
  signor_data <- read.table(signor,as.is = T,header = T,sep = '\t',fill = T,quote = "")
  signor_data_filtere_raw <- signor_data[signor_data$TypeA=='PROTEIN' & signor_data$TypeB=='PROTEIN',c('EntityA','IdA','EntityB','IdB','Effect','EffectMechanism','MechanismResidues')]
  signor_data_filtered <- signor_data_filtere_raw[tolower(signor_data_filtere_raw$EffectMechanism) %in% signor_interactions,]
  signor_table <- cbind(signor_data_filtered$EntityA,signor_data_filtered$IdA,signor_data_filtered$EntityB,signor_data_filtered$IdB,gsub(pattern = '[A-Za-z]+(\\d+)',replacement = '\\1',signor_data_filtered$MechanismResidues),tolower(signor_data_filtered$EffectMechanism),'SIGNOR')
  signor_table_unique <- unique(signor_table)
  colnames(signor_table_unique) <- c('GeneA','ProteinA','GeneB','ProteinB','Position','Type','Source')
  if(!include_auto){
    signor_table_unique <- signor_table_unique[signor_table_unique[,1] != signor_table_unique[,3] & signor_table_unique[,2] != signor_table_unique[,4],]
  }
  return(signor_table_unique)
}


#parse from Kinase_Substrate_Dataset
ks_parser <- function(ks,include_auto){
  ks_data <- read.csv(ks,header = T,as.is = T,skip = 2)
  ks_data_filtered <- ks_data[ks_data$SUB_ORGANISM == 'human' & ks_data$KIN_ORGANISM == 'human',]
  ks_data_table <- ks_data_filtered[,c('GENE','KIN_ACC_ID','SUB_GENE','SUB_ACC_ID','SUB_MOD_RSD')]
  ks_data_table$SUB_MOD_RSD <- gsub(pattern = '[A-Za-z]+',replacement = '',x = ks_data_table$SUB_MOD_RSD)
  ks_data_table_unique <- cbind(unique(ks_data_table),'phosphorylation','PhosphoSitePlus')
  colnames(ks_data_table_unique) <- c('GeneA','ProteinA','GeneB','ProteinB','Position','Type','Source')
  if(!include_auto){
    ks_data_table_unique <- ks_data_table_unique[ks_data_table_unique[,1] != ks_data_table_unique[,3] & ks_data_table_unique[,2] != ks_data_table_unique[,4],]
  }
  return(ks_data_table_unique)
}


#parse from biogrid
bg_parser <- function(bg,include_auto){
  bg_data <- read.table(bg,sep = '\t',fill = T,as.is = T,header = T,quote = '')
  bg_data_filtered <- bg_data[bg_data$Organism.Interactor.A=='9606' & bg_data$Organism.Interactor.B=='9606',c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
  bg_data_table_unique <- cbind(unique(bg_data_filtered),'PPI','BioGrid')
  colnames(bg_data_table_unique) <- c('GeneA','GeneB','Type','Source')
  if(!include_auto){
    bg_data_table_unique <- bg_data_table_unique[bg_data_table_unique[,1] != bg_data_table_unique[,2],]
  }
  return(bg_data_table_unique)
}

# parse from pathway commons
pw_parser <- function(pw,include_auto){
  pw_data <- read.table(pw,fill = T,as.is = T,header = F,sep = '\t')
  pw_data_table <- pw_data[pw_data$V2 %in% PC_interactions,]
  pw_data_table_unique <- unique(pw_data_table)
  if(!include_auto){
    pw_data_table_unique <- pw_data_table_unique[pw_data_table_unique$V1 != pw_data_table_unique$V3,]
  }
  pw_data_table_unique2 <- cbind(pw_data_table_unique$V1,pw_data_table_unique$V3,pw_data_table_unique$V2,'Pathway Commons')
  colnames(pw_data_table_unique2) <- c('GeneA','GeneB','Type','Source')
  return(pw_data_table_unique2)
}

expand_matrix <- function(m,ncol,sep = ';'){
  m[,ncol] <- as.character(m[,ncol])
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
interaction_table_site <- matrix(0,ncol=7,nrow = 0)
colnames(interaction_table_site) <- c('GeneA','ProteinA','GeneB','ProteinB','Position','Type','Source')
interaction_table_gene <- matrix(0,ncol=4,nrow = 0)
colnames(interaction_table_gene) <- c('GeneA','GeneB','Type','Source')

# SIGNOR
if(signor_file != ''){
  signor_table <- signor_parser(paste(interaction_data_path,signor_file,sep = '/'),signor_interactions = signor_interactions,include_auto = include_auto)
}else{
  signor_table <- NULL
}

#Kinase-substrate
if(ks_data_file != ''){
  ks_table <- ks_parser(paste(interaction_data_path,ks_data_file,sep = '/'),include_auto = include_auto)
}else{
  ks_table <- NULL
}

#BioGrid
if(bg_data_file != ''){
  bg_table <- bg_parser(paste(interaction_data_path,bg_data_file,sep = '/'),include_auto = include_auto)
}else{
  bg_table <- NULL
}

#Pathway Commons
if(pw_data_file != ''){
  pw_table <- pw_parser(paste(interaction_data_path,pw_data_file,sep = '/'),include_auto = include_auto)
}else{
  pw_table <- NULL
}

# KEGG(data already parsed from online database)
if(kegg_table_file != ''){
  kegg_table  <- read.csv(paste(output_path,kegg_table_file,sep = '/'),as.is = T)
  kegg_table2 <- unique(cbind(kegg_table[kegg_table$type %in% kegg_interactions, c(1,3,6)],'KEGG'))
  colnames(kegg_table2) <- c('GeneA','GeneB','Type','Source')
  if(!include_auto){
    kegg_table2 <- kegg_table2[kegg_table2$GeneA != kegg_table2$GeneB,]
  }
}
# Combine data from different sources
interaction_table_site <- unique(rbind(signor_table[signor_table[,'Position'] != '',],ks_table))
interaction_table_gene <- unique(rbind(signor_table[,c('GeneA','GeneB','Type','Source')],ks_table[,c('GeneA','GeneB','Type','Source')],bg_table,pw_table,kegg_table2))
interaction_table_site2 <- unique(expand_matrix(interaction_table_site,ncol = which(colnames(interaction_table_site)=='Position')))

interaction_table_gene2 <- apply(interaction_table_gene,2,as.character)
interaction_table_site3 <- apply(interaction_table_site2,2,as.character)
source_table <- rbind(cbind(interaction_table_site3[,6:7],paste(interaction_table_site3[,1],'_',interaction_table_site3[,2],'-',interaction_table_site3[,3],'_',interaction_table_site3[,4],'_',interaction_table_site3[,5],sep = '')),
                      cbind(interaction_table_gene2[,3:4],paste(interaction_table_gene2[,1],interaction_table_gene2[,2],sep = '-')))
colnames(source_table)[3] <- 'Implementation'

write.csv(unique(interaction_table_site3[,1:5]),paste(output_path,output_file_site,sep = '/'),row.names = F)
write.csv(unique(interaction_table_gene2[,1:2]),paste(output_path,output_file_gene,sep = '/'),row.names = F)
write.csv(source_table,paste(output_path,'source_table.csv',sep = '/'),row.names = F)




