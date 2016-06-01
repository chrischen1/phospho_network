library(glmnet)
library(methods)
library(randomForest)
library(parallel)
library(foreach)
library(doParallel)
library(igraph)

# Set parameters
result_path       = '~/Documents/workspace/phospho_network/example/script_files/output/test'
algos             = 'en'  # en for elastic net, rf for random forest
alphas            = seq(0,1) #required for elastic net
i_penalty         = T     # required for elastic net use different penalty based on heat diffusion?
ncore             = 2     # number of cores used
outerfold         = 5
innerfold         = 5
scale_method      = "0-1"   # 0-1 or "scale" 
directional       = T       # Used except pred_choose is flat or hf or all(hf would require undirectional version of hf_matrix. Should only upstream nodes in the pathway be considered?
pred_choose       = 'hf'    # method of choose different predictor : hf: by heat diffusion,connect:all connected nodes, direct: direct nodes, flat: all nodes in network, all: all data available by measurements
k                 = 0.001   # used if pred_choose is hf, a parameter to define the extend of predictors by keep nodes receive more heat than k*(heat_response_CNV)
max_level         = Inf     # used if pred_choose is up or connect, max level consdered for predictor selection

#LSF setting
cluster = F
queue   = 'short'
time    = '2:00'
mem     = '38000'

# Set inputs:
rna_filename          = '~/Documents/workspace/phospho_network/example/script_files/rna_processed.csv'
cnv_filename          = '~/Documents/workspace/phospho_network/example/script_files/cnv_processed.csv'
mut_filename          = '~/Documents/workspace/phospho_network/example/script_files/mutation_matrix_misc.csv'
mis_mut_filename      = '~/Documents/workspace/phospho_network/example/script_files/mutation_missense.csv'
heat_influence_file   = '~/Documents/workspace/phospho_network/example/script_files/network/heat_influence_test.csv' #used if pred_choose is hf
network_file          = '~/Documents/workspace/phospho_network/example/script_files/network/network_test.csv'

# target value input:
mdata_filename        = '~/Documents/workspace/phospho_network/example/script_files/rppa_processed.csv'

test_name='ARAF_299'

# args <- commandArgs(trailingOnly = TRUE)
# test_name <- gsub('+',';',args[1],fixed = T)
# config_filename <- args[2]
result_outfiles   = c(paste(test_name,'_','result_matrix.csv',sep = ''),paste(test_name,'_','beta_matrix.csv',sep = ''))
# source(config_filename)

# Functions:
# create vector for predictors with heat by heat diffusion
# test_node is in protein_site or gene_site format, both test_gene and test_node must in the constructed pathway network
predictor_heat <- function(test_node,hi_matrix,k = 0.01,test_gene = NULL){
  if(is.null(test_gene)){
    test_gene <- unique(gsub('(\\w+)_\\w+','\\1',test_node))
  }
  #create heat flow by heat influence and test_node
  heat_input <- matrix(0,nrow(hi_matrix),nrow(hi_matrix))
  colnames(heat_input) <- colnames(hi_matrix)
  rownames(heat_input) <- colnames(hi_matrix)
  # test_node <- test_node[test_node %in% names(diag(heat_input))]
  diag(heat_input)[test_node] <- 1/length(test_node)
  heat_flow <- hi_matrix %*% heat_input
  heat_flow2 <- apply(heat_flow,1,sum)
  predictors <- heat_flow2[heat_flow2 >= min(k*heat_flow2[test_gene])]
  return(predictors)
}

predictor_construct <- function(test_node_valid,pred_choose,network){
  if (pred_choose == 'connect'){
    bfs_dist <- get.all.shortest.paths(graph_from_adjacency_matrix(network),from = test_node_valid)$nrgeo
    predictors_name <- rownames(network)[bfs_dist > 0 & bfs_dist <= max_level]
  }else if (pred_choose == 'flat'){
    predictors_name <- rownames(network)
  }else if (pred_choose == 'direct'){
    test_gene <- strsplit(test_node_valid,'_')[[1]][1]
    interact_gene <- c(test_gene,unique(grep('_',rownames(network)[network[,test_gene]>0],value = T,invert = T)))
    predictors_name <- rownames(network)[unlist(lapply(strsplit(rownames(network),'_'),function(x)x[1])) %in% interact_gene]
  }
  predictors <- rep(1,length(predictors_name))
  names(predictors) <- predictors_name
  return(predictors_name)
}

# given a site_id, return a list with list$x are the predictors, list$y are true responses and list$penalty for individual penalty
data_prepare <- function(test_name,mdata,mut_data,cnv_data,rna_data,mis_mut,predictors,i_penalty,scale_method = '0-1'){
  # add rna info
  rna_predictors <- predictors[grep('_RNA',names(predictors))]
  rna_symbol <- gsub('(.+)_RNA','\\1',names(rna_predictors))
  rna_symbol <- rna_symbol[rna_symbol %in% rownames(rna_data)]
  rna_x <- t(rna_data[rna_symbol,])
  if(ncol(rna_x)>0){
    colnames(rna_x) <- paste(colnames(rna_x),'RNA',sep = '_')
  }
  # add cnv info
  cnv_predictors <- predictors[grep('_CNV',names(predictors))]
  cnv_symbol <- gsub('(.+)_CNV','\\1',names(cnv_predictors))
  cnv_symbol <- cnv_symbol[cnv_symbol %in% rownames(cnv_data)]
  cnv_x <- t(cnv_data[cnv_symbol,])
  if(ncol(cnv_x)>0){
    colnames(cnv_x) <- paste(colnames(cnv_x),'CNV',sep = '_')
  }
  # add mutation information
  mut_predictors <- predictors[grep('_mutation_misc',names(predictors))]
  mut_symbol <- gsub('(.+)_mutation_misc','\\1',names(mut_predictors))
  mut_symbol <- mut_symbol[mut_symbol %in% rownames(mut_data)]
  mut_x <- t(mut_data[mut_symbol,])
  valid_mut <- apply(mut_x,2,function(x)sum(x)>1)
  mut_x <- mut_x[,valid_mut,drop = F]
  if(sum(valid_mut) >= 1){
    colnames(mut_x) <- paste(colnames(mut_x),'mutation_misc',sep = '_')
  }
  mis_mut_predictors <- names(predictors)[names(predictors) %in% mis_mut$missense_names]
  mis_mut_x <- matrix(0,ncol = length(mis_mut_predictors),nrow = ncol(mdata))
  colnames(mis_mut_x) <- mis_mut_predictors
  rownames(mis_mut_x) <- colnames(mdata)
  if(length(mis_mut_predictors)>0){
    mis_mut_table <- mis_mut[mis_mut$missense_names %in% mis_mut_predictors,]
    for(i in 1:nrow(mis_mut_table)){
      mis_mut_x[gsub('-','.',mis_mut_table$sample[i]),mis_mut_table$missense_names[i]] <- 1
    }
  }
  intersect_samples <- Reduce(intersect,list(rownames(rna_x),rownames(cnv_x),rownames(mut_x),rownames(mis_mut_x)))
  data_x <- cbind(rna_x[intersect_samples,],cnv_x[intersect_samples,],mut_x[intersect_samples,],mis_mut_x[intersect_samples,])
  data_x <- apply(data_x,2,as.numeric)
  if(ncol(data_x) <= 2){return (NULL)}
  colnames(data_x) <- c(colnames(rna_x),colnames(cnv_x),colnames(mut_x),colnames(mis_mut_x))
  data_x_valid <- apply(data_x,2,function(x){x[is.na(x)] <- mean(x,na.rm=T);return(x)})
  pen_vec <- abs(log(predictors[colnames(data_x_valid)])+0.01)
  if(scale_method == '0-1'){
    data_x_norm <- apply(data_x_valid,2,data_norm)
    data_y_norm <- data_norm(as.numeric(t(mdata[test_name,intersect_samples])))
  }else{
    data_x_norm <- apply(data_x_valid,2,scale)
    data_y_norm <- scale(as.numeric(t(mdata[test_name,intersect_samples])))
  }
  # return x and y
  data_model <- list()
  data_model$x <- data_x_norm
  data_model$y <- data_y_norm
  if(max(pen_vec) == min(pen_vec) | !i_penalty){
    penalty <- predictors/predictors
  }else{
    penalty <- (pen_vec-min(pen_vec))/(max(pen_vec)-min(pen_vec))
  }
  data_model$penalty <- penalty
  data_model$site_id <- test_name
  return (data_model)
}

# prepare a list of sample for nested CV, list$train is a list of training indices and list$train is a list of testing indices
sample_prepare <- function(nsample,fold){
  test_list       <- list()
  train_list      <- list()
  sample_list     <- 1:nsample
  sample_list_tmp <- sample(sample_list,size = floor(length(sample_list)/outerfold)*outerfold)
  test_size <- floor(nsample/outerfold)
  for (i in 1:fold){
    test_list[[i]] <- sample_list_tmp[1:test_size]
    sample_list_tmp <- sample_list_tmp[-(1:test_size)]
    train_list[[i]] <- sample_list[-test_list[[i]]]
  }
  return(list('train' = train_list,'test' = test_list))
}

inner_rf <- function(model_data,train_set,innerfold){
  datax <- model_data$x[train_set,]
  datay <- model_data$y[train_set]
  # inner loop
  cv_model <- rfcv(datax,datay,cv.fold = innerfold)
  nvar = cv_model$n.var[which.min(cv_model$error.cv)]
  best_model <- randomForest(x = datax,y = datay,mtry = nvar,importance = T)
}

inner_en <- function(model_data,train_set,innerfold,alphas = alphas){
  datax <- model_data$x[train_set,]
  datay <- model_data$y[train_set]
  # inner loop
  best_alpha <- 0
  best_err   <- Inf
  # inner loop
  for(alpha in alphas){
    cv_model <- cv.glmnet(datax,datay,nfolds = innerfold,penalty.factor = model_data$penalty,alpha = alpha,grouped=FALSE)
    return_err <- min(cv_model$cvm)
    if(return_err < best_err){
      best_alpha <- alpha
      best_err   <- return_err
    }
  }
  cv_model <- cv.glmnet(datax,datay,nfolds = innerfold,penalty.factor = model_data$penalty,alpha = best_alpha,grouped=FALSE)
  best_model <- glmnet(datax,datay,alpha = best_alpha,lambda = cv_model$lambda.min,penalty = model_data$penalty)
  best_model$alpha <- best_alpha
  return(best_model)
}

inner_svm <- function(model_data,train_set,innerfold,alphas = alphas){
  kernels <- c('linear','polynomial','radial','sigmoid')
  datax <- model_data$x[train_set,]
  datay <- model_data$y[train_set]
  # inner loop
  best_model <- NULL
  best_err   <- Inf
  # inner loop
  for(k in kernels){
    cv_model <- svm(datax,datay,cross = innerfold,kernel = k)
    return_err <- cv_model$tot.MSE
    if(return_err < best_err){
      best_model <- cv_model
      best_err   <- return_err
    }
  }
  return(best_model)
}

# Nested cross validation: outer for performance evaluation, inner for best htperparameter search
nestcv <- function(model_data,ncore=2,outerfold,innerfold,algos,alphas = seq(0,1,0.05)){
  sample_list <- sample_prepare(nsample = nrow(model_data$x),fold = outerfold)
  result_matrix <- matrix(0,nrow = 0,ncol = 6)
  colnames(result_matrix) <- c('gene_site', 'test_set', 'parameter', 'predict_value', 'true_value', 'best_outer_q2')
  beta_matrix_raw   <- matrix(0,nrow = length(model_data$penalty),ncol = outerfold)
  colnames(beta_matrix_raw) <- c(paste('test',1:outerfold,sep = '_'))
  rownames(beta_matrix_raw) <- names(model_data$penalty)
  best_models <- list()
  train_list <- sample_list$train 
  test_list  <- sample_list$test
  # outer loop
  cl<-makeCluster(ncore)
  registerDoParallel(cl)
  if (algos == 'rf'){
    best_models <- foreach(i = 1:outerfold, 
                           .combine = list,.multicombine = TRUE,
                           .export =c('outerfold','inner_rf'),
                           .packages = 'randomForest')  %dopar% inner_rf(model_data,train_list[[i]],innerfold)
  }else if (algos == 'en'){
    best_models <- foreach(i = 1:outerfold, 
                           .combine = list,.multicombine = TRUE,
                           .export =c('outerfold','inner_en'),
                           .packages = 'glmnet')  %dopar% inner_en(model_data,train_list[[i]],innerfold,alphas)
  }
  stopCluster(cl)
  #end on outer loop
  for(i in 1:length(best_models)){
    test_x <- model_data$x[test_list[[i]],]
    if(outerfold == nrow(model_data$x)){
      test_x <- t(test_x)
    }
    outer_test_pred <- predict(best_models[[i]],test_x)
    #specific for algos
    if (algos == 'rf'){
      parameter <- best_models[[i]]$mtry
      beta_eval <- best_models[[i]]$importance[,1]
    }else if (algos == 'en'){
      parameter <- best_models[[i]]$alpha
      beta_eval <- best_models[[i]]$beta[,1]
    }
    result_matrix <- rbind(result_matrix,cbind(model_data$site_id,test_list[[i]],parameter,outer_test_pred,model_data$y[test_list[[i]]],0))
    beta_matrix_raw[names(beta_eval),i] <- beta_eval
  }
  pred_all <- as.numeric(result_matrix[,'true_value'])
  outer_q2 <- 1 - sum((as.numeric(result_matrix[,'predict_value'])-pred_all)^2)/sum((pred_all-mean(pred_all))^2)
  result_matrix[,'best_outer_q2'] <- rep(outer_q2,nrow(result_matrix))
  beta_matrix_raw2 <- beta_matrix_raw[apply(beta_matrix_raw,1,function(x)sum(abs(x))>0),]
  if(nrow(beta_matrix_raw2)>0){
    beta_matrix <- cbind('gene_site' = model_data$site_id,'predictor' = rownames(beta_matrix_raw2),beta_matrix_raw2)
  }
  rownames(beta_matrix) <- NULL
  return(list('matrix' = result_matrix,'beta_matrix' = beta_matrix))
}

data_norm <- function(x){
  if(max(x) != min(x)){return ((x-min(x))/(max(x)-min(x)))}
  else{return(rep(1,length(x)))}
}

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
mdata     <- read.csv(mdata_filename,row.names = 1)
rna_data  <- read.csv(rna_filename,as.is = T,row.names = 1)
mut_data  <- read.csv(mut_filename,as.is = T,row.names = 1)
cnv_data  <- read.csv(cnv_filename,as.is = T,row.names = 1)
mis_mut   <- read.csv(mis_mut_filename,as.is = T,row.names = 1)
# import network
network   <- as.matrix(read.csv(network_file,row.names = 1))
if(!directional){
  network <- get.adjacency(graph_from_adjacency_matrix(network,mode='undirected'),sparse=F)
}

#filter response avaliable in the network
node_all <- rownames(network)

# regression alalysis
result_matrix_all <- matrix(0,nrow = 0,ncol = 6)
beta_matrix_all   <- matrix(0,nrow = 0,ncol = 2+outerfold)
colnames(beta_matrix_all)   <- c('gene_site','predictor',paste('test',1:outerfold,sep = '_'))
colnames(result_matrix_all) <- c('gene_site', 'test_set', 'parameter', 'predict_value', 'true_value', 'best_outer_q2')

print(test_name)
response_list <- target_select(rownames(mdata))
test_node <- strsplit(response_list[response_list[,1] == test_name,2],split = ';')[[1]]
test_node_valid <- test_node[test_node %in% node_all]
if(length(test_node_valid) == 0){
  test_node_valid <- unique(gsub('(\\w+)_\\w+','\\1',test_node))
  test_node_valid <- test_node_valid[test_node_valid %in% node_all]
}
if(length(test_node_valid)>0){
  if(pred_choose == 'hf'){
    hi_matrix <- as.matrix(read.csv(heat_influence_file,row.names = 1))
    predictors <- predictor_heat(test_node_valid,hi_matrix,k=k)
  }else {
    predictors <- predictor_construct(test_node_valid,pred_choose = pred_choose,network = network)
  }
  model_data <- data_prepare(test_name = test_name,mdata = mdata,mut_data = mut_data,cnv_data = cnv_data,rna_data = rna_data,mis_mut = mis_mut,
                             predictors = predictors,i_penalty=i_penalty,scale_method = scale_method)
  if (!is.null(model_data)){
    model_result  <- nestcv(model_data = model_data,outerfold = outerfold,innerfold = innerfold,algos = algos,ncore = ncore)
    result_matrix_all <- model_result$matrix
    beta_matrix_all   <- model_result$beta_matrix
  }
  write.csv(result_matrix_all,paste(result_path,result_outfiles[1],sep = '/'),row.names = F)
  write.csv(beta_matrix_all,paste(result_path,result_outfiles[2],sep = '/'),row.names = F)
}