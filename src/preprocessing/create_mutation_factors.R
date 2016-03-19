# This script take mutation info of different cell lines from COSMIC and CCLE database, and output 2 tables: 
# First is a table of mutation types: 0:no mutaton; 1: mutations in UTR or intron; 2:mutations on CDS
# Second one would be the mutation from CCLE dataset.

#parameters
script_path     = '~/Documents/workspace/phospho_network/script_files'
path_cosmic     = '~/Documents/workspace/phospho_network/RAWDATA/mutation_data/cosmic'
file_ccle       = '~/Documents/workspace/phospho_network/RAWDATA/mutation_data/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.csv'
mcf_file        = '~/Documents/workspace/phospho_network/RAWDATA/mutation_data/mutations_mcf10a.csv'
msdata_file     = 'msdata_processed.csv'


cellines_ccle   = c('HCC38','BT549','MDAMB231','HCC1806','SKBR3','place_holder1','HCC70','BT-20','HS578T','MCF7')
# assign values for different types of mutation, the name of lists in list would be the value for its elements
mutation_type_ccle = list('1'=c("3'UTR","Intron","5'UTR","5'Flank"),
                          '2'=c("Missense_Mutation","Frame_Shift_Del","In_Frame_Del","Splice_Site_SNP",
                                "Nonsense_Mutation","Frame_Shift_Ins","Silent","Splice_Site_Ins","De_novo_Start_OutOfFrame",
                                "Nonstop_Mutation","Splice_Site_DNP","In_Frame_Ins","Splice_Site_Del","De_novo_Start_InFrame",
                                "Stop_Codon_DNP","Stop_Codon_Ins","Start_Codon_Del"))
cellines_cosmic = c('HCC38','BT-549','MDA-MB-231','HCC1806','place_holder1','place_holder2','HCC70','BT-20','Hs-578-T','MCF7')
# we exclude the non-coding and synonymous data in cosmic, so the only type should be 2
mutation_type_cosmic = 2
mutation_type_add    = 2
out_filename         = 'mutation_matrix.csv'

# Internal parameter:
cell_col = 11:20

# Main body
# Construct mutation matrix based on the data from CCLE
msdata <- read.csv(paste(script_path,msdata_file,sep = '/'),as.is = T,header = T)
ms_proteins <- unique(msdata$gene_symbol)
ms_cellline <- colnames(msdata)[cell_col]
mutation_matrix <- matrix(0,ncol = length(cell_col),nrow = length(ms_proteins))
colnames(mutation_matrix) <- colnames(msdata)[cell_col]
rownames(mutation_matrix) <- ms_proteins

data_ccle_raw <- read.csv(file = file_ccle,header = T,fill = T,as.is = T)
ccle_cell_id <- as.character(lapply(strsplit(data_ccle_raw$Tumor_Sample_Barcode,split = '_'),function(x)x[[1]]))
ind_filtered_ccle <- ccle_cell_id %in% cellines_ccle & data_ccle_raw$Hugo_Symbol %in% ms_proteins
data_ccle <- data_ccle_raw[ind_filtered_ccle,]
ccle_cell_id_filtered <- ccle_cell_id[ind_filtered_ccle]

for (i in 1:nrow(data_ccle)){
  mutaton_type <- data_ccle$Variant_Classification[i]
  mutation_num <- names(mutation_type_ccle)[lapply(mutation_type_ccle,function(x)mutaton_type %in% x)==T]
  mutation_matrix[data_ccle$Hugo_Symbol[i],which(cellines_ccle == ccle_cell_id_filtered[i])] <- mutation_num
}

# add more information based on COSMIC data
cosmic_file_list <- list.files(path_cosmic)
for (file in cosmic_file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("data_cosmic")){
    data_cosmic <- read.csv(paste(path_cosmic,file,sep = '/'), header=TRUE, as.is = T)
    data_cosmic <- cbind(cell= strsplit(file,split = '.',fixed = T)[[1]][1],data_cosmic)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("data_cosmic")){
    temp_dataset <-read.csv(paste(path_cosmic,file,sep = '/'), header=TRUE, as.is = T)
    data_cosmic<-rbind(data_cosmic, cbind(cell = strsplit(file,split = '.',fixed = T)[[1]][1], temp_dataset))
    rm(temp_dataset)
  }
}
data_cosmic$cell <- as.character(data_cosmic$cell)
data_cosmic$Gene <- as.character(lapply(strsplit(data_cosmic$Gene,split = '_'),function(x)x[1]))
ind_filtered_cosmic <- data_cosmic$Gene %in% msdata$gene_symbol
data_cosmic_filtered <- data_cosmic[ind_filtered_cosmic,]

for (i in 1:nrow(data_cosmic_filtered)){
  mutation_matrix[data_cosmic_filtered$Gene[i],data_cosmic_filtered$cell[i]] <- mutation_type_cosmic
}

# add mutation data for MCF10A
data_add <- read.csv(mcf_file,header = F,as.is = T)
gene_names_add <- data_add$V1
gene_names_add_filtered <- gene_names_add[gene_names_add %in% ms_proteins]
mutation_matrix[gene_names_add_filtered,"MCF10A"] = mutation_type_add

# output result
write.csv(mutation_matrix,file = paste(script_path,out_filename,sep = '/'))

