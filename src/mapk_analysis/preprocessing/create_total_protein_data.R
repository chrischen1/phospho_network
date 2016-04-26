library(preprocessCore)

# parsing data of total protein
total_protein_path     = '/Users/xxx/Documents/workspace/phospho_network/processed_data'
script_dir             = '/Users/xxx/Documents/workspace/phospho_network/script_files'

total_protein_list     = c('RepAwithPhos_0.6IS_300sn_031215.csv','RepBwithPhos_0.6IS_300sn_031215.csv')
total_protein_outfile  = 'total_protein_processed.csv'
total_cell_col = 8:17

# helper funcitons:
# quantile normalize and average the colnumns given by cell_col
qnorm_avg <- function(msdata_raw,avg_col = cell_col){
  m <- msdata_raw[1,]
  data_qnorm <- t(normalize.quantiles(t(msdata_raw[,avg_col])))
  m[avg_col] <- apply(data_qnorm,2,mean)
  return(m)
}

# main body
id_list <- list()
total_protein_raw <- list()
total_protein_intersect <- list()
for (i in 1:length(total_protein_list)){
  data_raw <- read.csv(paste(total_protein_path,total_protein_list[i],sep = '/'),as.is = T, header = T)
  acc_id <- gsub(pattern = '\\w+\\|(.+)\\|\\w+',replacement = '\\1',data_raw$Protein.Id)
  total_protein_raw[[i]] <- cbind.data.frame(acc_id,data_raw)
  id_list[[i]] <- acc_id
}
intersect_ids <- Reduce(intersect, id_list)
for (i in 1:length(total_protein_raw)){
  total_protein <- total_protein_raw[[i]]
  protein_intersect <- total_protein[total_protein$acc_id %in% intersect_ids,]
  total_protein_intersect[[i]] <- protein_intersect[order(protein_intersect$acc_id),]
}
total_protein_processed <- total_protein_intersect[[1]]
for (i in 1:length(intersect_ids)){
  m_raw <- do.call(rbind,lapply(total_protein_intersect,function(x)x[i,]))
  m_avg <- qnorm_avg(m_raw,avg_col = total_cell_col)
  total_protein_processed[i,] <- m_avg
}

merge_id <- gsub('(\\w+)-\\d+','\\1',total_protein_processed$acc_id)
total_protein_processed2 <- cbind(merge_id,total_protein_processed)
write.csv(total_protein_processed2,file = paste(script_dir,total_protein_outfile,sep = '/'),row.names = F)
