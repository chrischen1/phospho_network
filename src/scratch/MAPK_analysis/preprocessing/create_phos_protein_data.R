# This Rscript would take a list of files(each corresponding to a sheet in excel file) as input, and do the following steps:
# 1.Remove any data with max score less than score_cutoff
# 2.Remove any data with default_num_quant less than quant_cutoff
# 3.Merge duplicated information(same uniprot ID, same position) by QNORM average
# 4.(Optional) Keep only intersect data between replicates A and B, merge them by QNORM average
# 5.Merge phospolation sites with correlations larger than cor_cutoff percent of correlations between replicates
# Finally, processed data would ouput as a csv file which name given by outfile_name

library(preprocessCore)
#Set parameters:
infile_dir      = "/Users/xxx/Documents/workspace/phospho_network/processed_data"
infile_list     = list(c('ReplicateA_pSTY_Summary_031315_1.csv','ReplicateA_pSTY_Summary_031315_2.csv',
                        'ReplicateA_pSTY_Summary_031315_3.csv','ReplicateA_pSTY_Summary_031315_4.csv'),
                       c('ReplicateB_pSTY_Summary_031315_1.csv','ReplicateB_pSTY_Summary_031315_2.csv',
                        'ReplicateB_pSTY_Summary_031315_3.csv','ReplicateB_pSTY_Summary_031315_4.csv'))
script_dir      = "/Users/xxx/Documents/workspace/phospho_network/script_files"
outfile_name    = 'msdata_processed.csv'
score_cutoff    = 19        #Remove any data with max score less than score_cutoff
quant_cutoff    = 1         #Remove any data with default_num_quant less than quant_cutoff
cor_cutoff      = 0.5       #Merge phospolation sites with correlations larger than cor_cutoff percent of correlations between replicates
filter_rep      = T     #Keep only intersect data between replicates A and B?

# Internal parameters
cell_col        = 10:19


# Functions:
# Integrate MS data into one table, keeping all information
integrate_msdata <- function(file_list,infile_dir){
  data1 <- read.csv(paste(infile_dir,file_list[1],sep = '/'),header = T,skip = 1,as.is = T)
  data3 <- read.csv(paste(infile_dir,file_list[3],sep = '/'),header = T,skip = 1,as.is = T)
  data_single <- rbind(data1,data3)
  data_single2 <- data_single[,-(14:16)]
  data2 <- read.csv(paste(infile_dir,file_list[2],sep = '/'),header = T,skip = 1,as.is = T)
  data4 <- read.csv(paste(infile_dir,file_list[4],sep = '/'),header = T,skip = 1,as.is = T)
  data_mult <- rbind(data2,data4)
  data_mult2 <- data_mult[,c(1,3:6,8,7,9:25)]
  colnames(data_mult2) <- colnames(data_single2)
  data_all <- rbind(data_single2,data_mult2)
  cell_names <- colnames(read.csv(paste(infile_dir,file_list[1],sep = '/'),nrow=1))[-(1:17)]
  colnames(data_all)[-(1:14)] = cell_names
  acc_id <- gsub(pattern = '\\w+\\|(.+)\\|\\w+',replacement = '\\1',data_all$Protein.Id)
  return(data_all[,c(2:7,13:24)])
}

# Extract and add uniprot ID
add_acc_id <- function(msdata_raw){
  acc_id <- gsub(pattern = '\\w+\\|(.+)\\|\\w+',replacement = '\\1',msdata_raw$Protein.Id)
  return(cbind('acc_id' = acc_id,msdata_raw))
}

# Remove any data with quanlity score less than score_cutoff(mutiple sites using the min among them)
filter_score <- function(msdata_raw, cutoff = score_cutoff){
  scores_list <- strsplit(msdata_raw$Max.Score,';')
  min_scores <- rep(0,length(scores_list))
  for (i in 1:length(scores_list)){
    min_scores[i] <- min(as.numeric(scores_list[[i]]))
  }
  msdata_raw$Max.Score <- min_scores
  return(msdata_raw[min_scores >= score_cutoff,])
}

# Remove any data with quant less than quant_cutoff
filter_quant <- function(msdata_raw, cutoff = quant_cutoff){
  return(msdata_raw[msdata_raw$default_num_quant >= quant_cutoff,])
}

#remove duplicates(same id, same site) by keeping QNORM avg of them, data averged are col in datacol
avg_dup <- function(msdata_raw,avg_col = cell_col){
  ids <- paste(msdata_raw$acc_id,msdata_raw$Site.Position,sep = '_')
  unique_id <- unique(ids)
  dup_index <- duplicated(ids)
  dup_id <- ids[dup_index]
  no_dup_index <- !(ids %in% dup_id)
  msdata_unique <- msdata_raw[no_dup_index,]
  for (i in 1:length(dup_id)){
    dup_m <- msdata_raw[ids == dup_id[i],]
    m <- qnorm_avg(dup_m,avg_col = cell_col)
    m$Max.Score <- min(dup_m$Max.Score)
    m$default_num_quant <- min(dup_m$default_num_quant)
    msdata_unique <- rbind(msdata_unique,m)
  }
  return(as.data.frame(msdata_unique))
}

# Merge data of replicates with QNORM average, input is a list of lists
merge_replicate <- function(data_raw_list,intersect_id_list,avg_col = cell_col){
  num_id <- length(intersect_id_list[[1]])
  num_rep <- length(data_raw_list)
  msdata_merge <- data_raw_list[[1]][num_id,]
  for (i in 1:num_id){
    #Average ith protein in uniqueID_list
    rep_m <- data_raw_list[[1]][1:num_rep,]
    for (j in 1:num_rep){
      rep_ind <- intersect_id_list[[j]][i]
      rep_m[j,] <- data_raw_list[[j]][rep_ind,]
    }
    m <- qnorm_avg(rep_m,avg_col = avg_col)
    m$Max.Score <- min(rep_m$Max.Score)
    m$default_num_quant <- min(rep_m$default_num_quant)
    msdata_merge[i,] <- m
  }
  return(as.data.frame(msdata_merge))
}

# Merge phospolation sites with correlations larger than cor_cutoff percent of correlations between replicates
merge_phos_sites <- function(data_merge,cutoff = cor_cutoff_value, cell_col = cell_col){
  all_proteins <- as.character(data_merge$acc_id)
  multi_site_protein <- all_proteins[duplicated(all_proteins)]
  data_merge_phos <- data_merge[!all_proteins %in% multi_site_protein,]
  unique_multi_site_protein <- unique(multi_site_protein)
  for (i in 1:length(unique_multi_site_protein)){
    protein_name <- unique_multi_site_protein[i]
    m <- data_merge[all_proteins == protein_name,]
    m_dist <- as.dist(1-cor(t(m[,cell_col]),method = 'spearman'))
    clust <- hclust(m_dist)
    cut <- cutree(clust,h = 1-cor_cutoff_value)
    for (j in 1:max(cut)){
      m_group <- m[cut==j,]
      m_group_qnorm <- qnorm_avg(m_group,avg_col = cell_col)
      data_merge_phos <- rbind(data_merge_phos,m_group_qnorm)
    }
  }
  return(data_merge_phos)
}

# geneate a list of lists, list[[i]] are the indices of intersect elements in ith replicate
# list[[i]][k] and list[[j]][k] would be the same unique ID
intersect_indices <- function(data_raw_list){
  id_list <- list()
  intersect_indices <- list()
  for (i in 1:length(data_raw_list)){
    id_list[[i]] <- paste(data_raw_list[[i]]$acc_id,data_raw_list[[i]]$Site.Position,sep = '_')
  }
  common_ids <- Reduce(intersect, id_list)
  for (i in 1:length(data_raw_list)){
    intersect_id_list <- c()
    for (j in 1:length(common_ids)){
      intersect_id_list[j] <- which(id_list[[i]] == common_ids[j])
    }
    intersect_indices[[i]] <- intersect_id_list
  }
  return (intersect_indices)
}

# helper funcitons:
# quantile normalize and average the colnumns given by cell_col
qnorm_avg <- function(msdata_raw,avg_col = cell_col){
  m <- msdata_raw[1,]
  data_qnorm <- t(normalize.quantiles(t(msdata_raw[,avg_col])))
  m[avg_col] <- apply(data_qnorm,2,mean)
  return(m)
}

# Calculate all correlations between replicates pairs(same protein, same position)
# If there are more than 2 replicates, low triangle of the correlation matrix would be taken in the whole correlations
replicate_cor <- function(data_raw_list,intersect_id_list,cell_col = cell_col){
  replicate_cor_value <- c()
  num_id <- length(intersect_id_list[[1]])
  num_rep <- length(data_raw_list)
  for(i in 1:num_id){
    m <- t(rbind(data_raw_list[[1]][intersect_id_list[[1]][i],cell_col],data_raw_list[[2]][intersect_id_list[[2]][i],cell_col]))
    cor_m <- cor(m)
    cors <- cor_m[lower.tri(cor_m)]
    replicate_cor_value <- c(replicate_cor_value,cors)
  }
  return(replicate_cor_value)
}

# Main body:
data_raw_list <- list()
for (i in 1:length(infile_list)){
  msdata_raw <- integrate_msdata(infile_list[[i]],infile_dir = infile_dir)
  msdata_raw <- add_acc_id(msdata_raw)
  msdata_raw_score <- filter_score(msdata_raw)
  msdata_raw_score_quant <- filter_quant(msdata_raw_score)
  msdata_raw_score_quant_nodup <- avg_dup(msdata_raw_score_quant)
  data_raw_list[[i]] <- msdata_raw_score_quant_nodup
}

if(length(infile_list) > 1){
  intersect_id_list <- intersect_indices(data_raw_list)
  data_merge <- merge_replicate(data_raw_list,intersect_id_list = intersect_id_list, avg_col = cell_col)
  cor_cutoff_value <- quantile(replicate_cor(data_raw_list,intersect_id_list = intersect_id_list,cell_col = cell_col),cor_cutoff)
  data_processed <- merge_phos_sites(data_merge, cutoff = cor_cutoff_value, cell_col = cell_col)
  if(!filter_rep){
    for (i in 1:length(infile_list)){
      data_processed <- rbind(data_processed,data_raw_list[[i]][-intersect_id_list[[i]],])
    }
  }
}else{
  data_processed <- data_raw_list[[1]]
}
merge_id <- gsub('(\\w+)-\\d+','\\1',data_processed$acc_id)
data_output <- cbind(merge_id,data_processed)
write.csv(data_output,paste(script_dir,outfile_name,sep = '/'),row.names = F)
