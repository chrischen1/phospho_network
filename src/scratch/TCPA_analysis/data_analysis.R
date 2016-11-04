# This Rscript would take a processed measured data file as input, and do the following steps:
# 1.Import the MS data and network file(as edge list) from csv files
# 2.Filter the data by keeping only proteins in the network 
# 3.Nodes 1 level further upstream could also be treated as predictors, but a penalty factor would be imposed
# 4.Construct the LASSO model: PeptideA_site1 ~ PeptideA_site1 + PeptideA_site2...
# 5.Evaluate the model with r-square and q-square
# input files:
# 1. measure_filename tot_prt_filename are processed quantitive data for phos or total protein
# 2. interaction_gene_file: info for gene interactions or gene interactions to build the network
# 3. interaction_site_file: addtional info for protein/site interactions or gene interactions
# 4. mutation_filename is info about mutations on each sample
# 5. target_gene_file: genes of interest used as response
# output files:
# 1.result_matrix is the prediction versus true value from Nested cross validation, it also contains the value for hyperparameters trained from inner loop
# 2.beta_matrix is the table for coefficient from different testsets.

library(glmnet)
library(methods)
library(igraph)

# Set parameters
script_dir        = '~/Documents/workspace/phospho_network/script_files'
results_dir       = '~/Documents/workspace/phospho_network/script_files/analysis_rppa/mapk'
network_dir       = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed'

measure_filename      = '~/Documents/workspace/phospho_network/script_files/msdata_processed.csv'
tot_prt_filename      = '~/Documents/workspace/phospho_network/script_files/total_protein_processed.csv'
interaction_site_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/interaction_table_site_all.csv'
interaction_gene_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/network_all_gene.csv'
mutation_filename     = '~/Documents/workspace/phospho_network/script_files/mutation_matrix.csv'
target_gene_file      = '~/Documents/workspace/phospho_network/script_files/analysis_rppa/response_list.txt'

result_outfiles   = c('result_matrix.csv','beta_matrix.csv')

data_col          = 11:20
total_data_col    = 9:18
penalty_site      = 0.01
penalty_prot      = 1
penalty_gene      = 2
penalty_noevd     = 3
penalty_tot       = 1
penalty_mut       = 1

alphas            = seq(0,1,0.05)
lambdas           = c()
indirect_predictor= 0
num_test          = 1 #number of test sets for outer loop
outerfold         = 10
innerfold         = 9

# Functions:
# given a site_id, return a list with list$x are the predictors, and list$y are true responses
data_prepare <- function(ind,mdata,totdata = NULL,network,interaction_site= NULL,interaction_prot = NULL,interaction_gene= NULL,data_col = data_col,
                         total_data_col = total_data_col,mutation_matrix = NULL,penalty_site= 0.01,penalty_prot=1,penalty_gene=2, penalty_noevd=3,penalty_tot=1,penalty_mut = 1){
  data_model <- list(x=c(),y=c())
  gene_name <- mdata$gene_symbol[ind]
  prot_name <- mdata$merge_id[ind]
  site_list <- as.character(as.numeric(strsplit(mdata$Site.Position[ind],split = ';')[[1]]))
  all_predictors_genes <- rownames(network)[network[,gene_name] > 0]
  all_predictors_data  <- mdata[mdata$gene_symbol %in% all_predictors_genes,]
  penalty_vector <- rep(penalty_noevd,nrow(all_predictors_data))
  if(nrow(all_predictors_data) < 1){
    return(NULL)
  }else{
    if(!is.null(interaction_gene)){
      gene_evidence <- unique(interaction_gene[interaction_gene$geneB == gene_name,'geneA'])
      penalty_vector[which(all_predictors_data$gene_symbol %in% gene_evidence)] <- penalty_gene
    }
    if (!is.null(interaction_prot)){
      prot_evidence <- unique(interaction_prot[interaction_prot$protB == prot_name,'protA'])
      penalty_vector[which(all_predictors_data$merge_id %in% prot_evidence)] <- penalty_prot
    }
    if (!is.null(interaction_site)){
      site_evidence <- interaction_site[interaction_site$geneB == gene_name,]
      for (i in nrow(all_predictors_data)){
        p_gene_name <- all_predictors_data$gene_symbol[i]
        if (p_gene_name %in% site_evidence$geneA){
          p_site_list <- as.character(as.numeric(strsplit(all_predictors_data$Site.Position[i],';')[[1]]))
          if(length(intersect(site_evidence$site[site_evidence$geneA == p_gene_name],p_site_list))>1){
            penalty_vector[i] <- penalty_site
          }
        }
      }
    }
    x <- t(all_predictors_data[,data_col])
    colnames(x) <- paste(all_predictors_data$acc_id,all_predictors_data$Site.Position,sep = '_')
    if(is.null(totdata)){
      if(gene_name %in% totdata$gene_symbol){
        tot_protein <- as.numeric(totdata[totdata[,'gene_symbol'] == gene_name,total_data_col])
        x <- cbind(x,'tot_protein' = tot_protein)
        penalty_vector <- c(penalty_vector,penalty_tot)
      }
    }
    if(!is.null(mutation_matrix)){
      if(gene_name %in% colnames(mutation_matrix)){
        x <- cbind(x,'mutation' = as.numeric(mutation_matrix[gene_name,]))
        penalty_vector <- c(penalty_vector,1)
      }
    }
    x2 <- apply(x,2,as.numeric)
    if(ncol(x2) < 2){return (NULL)}
    rownames(x2) <- rownames(x)
    colnames(x2) <- colnames(x)
    data_model$x <- x2
    data_model$y <- as.numeric(mdata[ind,data_col])
    data_model$penalty <- penalty_vector
    data_model$site_id <- paste(gene_name,mdata$Site.Position[ind],sep = '_')
  }
  return (data_model)
}

# Nested cross validation: alpha and lambda are vectors for grid search, if lambda is set to NULL, it would train alpha only and let glmnet to decide lambda
# inner loop would use a leave-1-out CV, outer loop would have ntest for test set.
nestcv <- function(model_data,alphas,lambdas = NULL,outerfold = 3,innerfold = 10){
  penalty         <- model_data$penalty
  test_list       <- list()
  sample_list     <- 1:nrow(model_data$x)
  sample_list_tmp <- sample(sample_list,size = floor(length(sample_list)/outerfold)*outerfold)
  for (i in 1:outerfold){
    begins <- floor(length(sample_list)/outerfold)*i
    ends   <- floor(length(sample_list)/outerfold)*(i+1)-1
    test_list[[i]] <- sample_list_tmp[begins:ends]
  }
  result_matrix <- matrix(0,nrow = 0,ncol = 6)
  colnames(result_matrix) <- c('gene_site','test_set','alpha','predict_value','true_value','best_outer_q2')
  beta_matrix_raw <- matrix(0,nrow = 0,ncol = ncol(model_data$x))
  outer_pred_list <- rep(0,nrow(model_data$x))
  # outer loop
  for (i in 1:outerfold){
    outer_train_set <- sample_list[-test_list[[i]]]
    outer_test_set  <- test_list[[i]]
    best_alpha <- 0
    best_err   <- Inf
    # inner loop
    for(alpha in alphas){
      cv_model <- cv.glmnet(model_data$x[-test_list[[i]],],model_data$y[-test_list[[i]]],nfolds = innerfold,penalty.factor = model_data$penalty,alpha = alpha,grouped=FALSE)
      return_err <- min(cv_model$cvm)
      if(return_err < best_err){
        best_alpha <- alpha
        best_model <- cv_model$glmnet.fit
      }
    }
    #end on inner loop
    # use best alpha on test data
    test_x <- model_data$x[outer_test_set,]
    if(outerfold == nrow(model_data$x)){
      test_x <- t(test_x)
    }
    outer_test_pred <- predict(best_model,newx = test_x,s = cv_model$lambda.min)
    result_matrix <- rbind(result_matrix,cbind(rep(model_data$site_id,length(outer_test_set)),outer_test_set,rep(best_alpha,length(outer_test_set)),outer_test_pred,model_data$y[outer_test_set],rep(0,length(outer_test_set))))
    
    beta_matrix_raw <- rbind(beta_matrix_raw,best_model$beta[,which.max(best_model$dev.ratio)])
  }
  pred_all <- as.numeric(result_matrix[,'true_value'])
  outer_q2 <- 1 - sum((as.numeric(result_matrix[,'predict_value'])-pred_all)^2)/sum((pred_all-mean(pred_all))^2)
  
  beta_matrix <- t(rbind(rep(model_data$site_id,ncol(beta_matrix_raw)),colnames(model_data$x),beta_matrix_raw))
  colnames(beta_matrix) <- c('gene_site','predictor',paste('test',1:outerfold,sep = '_'))
  result_matrix[,'best_outer_q2'] <- rep(outer_q2,nrow(result_matrix))
  return(list('matrix' = result_matrix,'best_outer_q2' = outer_q2,'beta_matrix' = beta_matrix))
}

# return a matrix with all upstream nodes that are one level further(all_level=F) or every possible upstream nodes(all_level=F)
up_network <- function(network,all_level=F,penalty = 2){
  network_m <- apply(network,2,as.numeric)
  if(!all_level){
    upstream <- network_m %*% network_m
    upstream[upstream >= 1] <- penalty
    upstream[network == 1] <- 1
  }else{
    upstream <- network
    level_up <- network_m %*% network_m
    current_penalty <- penalty
    while(sum(level_up) > 0){
      upstream[level_up >= 1 & upstream == 0] <- current_penalty
      current_penalty <- current_penalty*penalty
      level_up <- level_up %*% level_up
    }
  }
  rownames(upstream) <- colnames(upstream)
  return(upstream)
}

cal_q2 <- function(true,pred){
  pred <- as.numeric(pred)
  true <- as.numeric(true)
  return(1 - sum((pred-true)^2)/sum((true-mean(true))^2))
}

expand_matrix <- function(m,ncol,sep = '; '){
  m2 <- m[0,]
  collision <- c()
  for(i in 1:nrow(m)){
    var_split <- strsplit(m[i,ncol],split = sep)[[1]]
    if(length(var_split) <= 1){
      m2 <- rbind(m2,m[i,])
      collision <- c(collision,'')
    }else if(length(var_split) > 1){
      # avoid 'perfect prediction'
      new_table <- matrix(rep(as.character(m[i,]),length(var_split)),nrow = length(var_split),byrow = T)
      new_table[,ncol] <- var_split
      colnames(new_table) <- colnames(m2)
      m2 <- rbind(m2,new_table)
      collision <- c(collision,rep(m[i,ncol],length(var_split)))
    }
  }
  m3 <- cbind(m2,collision)
  m3[,'collision'] <- as.character(m3[,'collision'])
  return(m3)
}

# Main
target_list <- unlist(read.table(target_gene_file,as.is = T))
tar_genes  <- gsub('(\\w+)_.+','\\1',target_list)
mutation_matrix <- read.csv(mutation_filename,as.is = T,row.names = 1)
mdata_raw <- read.csv(measure_filename,as.is = T,header = T)
totdata  <- read.csv(tot_prt_filename,as.is = T)
interaction_site_raw <- read.csv(interaction_site_file,as.is = T)
interaction_gene_raw <- read.csv(interaction_gene_file,as.is = T)
if(sum(grep(';',mdata_raw$gene_symbol)) > 0){
  mdata_raw <- expand_matrix(mdata_raw,ncol = which(colnames(mdata_raw) == 'gene_symbol'))
}
all_genes <- intersect(unique(mdata_raw$gene_symbol),unique(unlist(interaction_gene_raw)))
mdata <- mdata_raw[mdata_raw$gene_symbol %in% all_genes,]
rownames(mdata) <- 1:nrow(mdata)

interaction_site <- interaction_site_raw[interaction_site_raw$geneA %in% all_genes & interaction_site_raw$geneB %in% all_genes & interaction_site_raw$geneA != interaction_site_raw$geneB,]
interaction_gene <- interaction_gene_raw[interaction_gene_raw$geneA %in% all_genes & interaction_gene_raw$geneB %in% all_genes & interaction_gene_raw$geneA != interaction_gene_raw$geneB,]

g <- graph_from_edgelist(as.matrix(interaction_gene))
network <- as.matrix(as_adjacency_matrix(g))[,tar_genes]

result_matrix <- matrix(0,nrow = 0,ncol = 6)
beta_matrix   <- matrix(0,nrow = 0,ncol = 2+outerfold)
colnames(beta_matrix)   <- c('gene_site','predictor',paste('test',1:outerfold,sep = '_'))
colnames(result_matrix) <- c('gene_site', 'test_set', 'alpha', 'predict_value', 'true_value', 'best_outer_q2')

# select target index:
gene_site <- paste(mdata$gene_symbol,mdata$Site.Position,sep = '_')
tar_index <- which(gene_site %in% target_list)
if(length(gene_site[tar_index]) > length(unique(gene_site[tar_index]))){
  warning('duplicate unique id found in target')
}

# iterate of all phos sites on the dataset
for (ind in tar_index){
  model_data  <- data_prepare(ind = ind,mdata = mdata,totdata = totdata,network = network,interaction_site = interaction_site, mutation_matrix = mutation_matrix,
                              data_col = data_col,total_data_col = total_data_col,penalty_site = penalty_site,penalty_prot = penalty_prot,penalty_gene = penalty_gene,penalty_noevd = penalty_noevd,penalty_tot = penalty_tot,penalty_mut = penalty_mut)
  if (!is.null(model_data)){
    model_result  <- nestcv(model_data,alphas = alphas,outerfold = outerfold,innerfold = innerfold)
    result_matrix <- rbind(result_matrix,model_result$matrix)
    beta_matrix   <- rbind(beta_matrix,model_result$beta_matrix)
  }
}

write.csv(result_matrix,paste(results_dir,result_outfiles[1],sep = '/'),row.names = F)
write.csv(beta_matrix,paste(results_dir,result_outfiles[2],sep = '/'),row.names = F)
