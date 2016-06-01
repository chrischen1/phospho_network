# usage: regression_mcore_parallel_control.R config_file_name 
# This Rscript would take a processed measured data file as input, and do the following steps:
# 1.Import the protein data and network file(as edge list) from csv files
# 2.Import matrix for RNA, CNV and mutation as predictors
# 3.Select predictors based on heat diffusion from response in the network
# 4.Construct the LASSO model: PeptideA_site1 ~ PeptideA_site1 + PeptideA_site2...
# 5.Evaluate the model with r-square and q-square
# input files:
# 1. measure_filename are processed quantitive data for phos protein
# 2. heat_influence_file: file for heat influence matrix
# 3. network_file: file for pathway network
# 4. mut_filename/rna_filename/cnv_filename is info about mutations/expression/copy_number info on each sample
# output files:
# 1.result_matrix is the prediction versus true value from Nested cross validation, it also contains the value for hyperparameters trained from inner loop
# 2.beta_matrix is the table for coefficient(elastic net) or importance(random forest) from different testsets.

# the order of command line args given to regression_mcore_parallel_regression.R:
# test_name,rna_filename,cnv_filename,mut_filename,mis_mut_filename,heat_influence_file,network_file,mdata_filename,result_path
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
                 'Rscript regression.R',test_name,config_filename))
  }else{
    system(paste('Rscript regression.R',test_name,config_filename))
  }
}