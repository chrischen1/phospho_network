result_path = '/n/scratch2/cc400/align_star/bam_2pass'
output_file = '~/Documents/workspace/phospho_network/RNAseq/processed_data/star/star_log.txt'

log_files <- grep('Log.final.out',list.files(result_path),value = T)
log_table    <- NULL
for (i in log_files){
  l <- readLines(paste(result_path,i,sep = '/'))
  l_list <- strsplit(l,'|\t',fixed = T)
  l_val <- sapply(l_list,function(x)x[2])
  names(l_val) <- gsub(' +',' ',sapply(l_list,function(x)x[1]))
  log_table <- cbind(log_table,l_val[!is.na(l_val)][-(1:4)])
}
colnames(log_table) <- gsub('_qcLog.final.out','\\1',log_files)
write.table(log_table,output_file,sep = '\t',row.names = T)