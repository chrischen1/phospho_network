args <- commandArgs(trailingOnly = TRUE)
config_filename <- args[1]

# Functions:
target_select <- function(mdata_names){
  mdata_p <- grep('\\w+_\\w+',mdata_names,value = T)
  mdata_plist <- strsplit(mdata_p,split = '_')
  mdata_pgenes <- strsplit(unlist(lapply(mdata_plist,function(x)x[1])),split = ';')
  mdata_psites <- strsplit(unlist(lapply(mdata_plist,function(x)x[2])),split = ';')
  mdata_pnodes <- c()
  for(i in 1:length(mdata_pgenes)){
    new_nodes=expand.grid(mdata_pgenes[[i]],mdata_psites[[i]])
    mdata_pnodes <- c(mdata_pnodes,paste(paste(new_nodes$Var1,new_nodes$Var2,sep = '_'),collapse = ';'))
  }
  return(cbind(mdata_p,mdata_pnodes))
}

# Main
# import dataset
source(config_filename)
mdata <- read.csv(mdata_filename,row.names = 1)
dir.create(file.path(result_path, ''), showWarnings = FALSE)

#filter response avaliable in the network
response_list <- target_select(rownames(mdata))
response_list <- gsub(';','+',response_list)

for(test_name in response_list[,1]){
  if(cluster){
    print(paste(test_name,'submitted'))
    system(paste('bsub -q',queue,'-W',time,'-R',paste('"rusage[mem=',mem,']"',sep = ''),
                 '-o',paste(result_path,'/',test_name,'_job.out',sep = ''),
                 '-e',paste(result_path,'/',test_name,'_job.err',sep = ''),
                 'Rscript predictor_analysis.R',test_name,config_filename))
  }else{
    system(paste('Rscript predictor_analysis.R',test_name,config_filename))
  }
}

