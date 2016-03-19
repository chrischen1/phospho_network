# This Rscript would take a processed MS data file(data_processing.R) as input, and do the following steps:
# 1.Import the MS data and network file(as an adjacency matrix) from csv files
# 2.Filter the data by keeping only proteins in the network 
# 3.Nodes 1 level further upstream could also be treated as predictors, but a penalty factor would be imposed
# 4.Construct the LASSO model: PeptideA_site1 ~ PeptideA_site1 + PeptideA_site2...
# 5.Evaluate the model with r-square and q-square
# input files:
# 1. msdata_filename tot_prt_filename are processed MS data for phos or total protein
# 2. network_filename: network for mapk pathway as adj matrix for literal including kegg
# 3. prot_table_file gene_table_file addtional info for protein/site interactions or gene interactions
# 4. mutation_filename is info about mutations on each cell lines

library(glmnet)

# Set parameters
script_dir        = '~/Documents/workspace/phospho_network/script_files'
results_dir       = '~/Documents/workspace/phospho_network/script_files/analysis_results/mapk_merge'
network_dir       = '~/Documents/workspace/phospho_network/script_files/mapk_analysis'


msdata_filename   = 'msdata_processed.csv'
tot_prt_filename  = 'total_protein_processed.csv'
network_filename  = 'mapk_network_literal.csv'
prot_table_file   = 'mapk_prot_table.csv'
gene_table_file   = 'mapk_gene_table.csv'
mutation_filename = 'mutation_matrix.csv'
alphas            = seq(0,1,0.05)
lambdas           = c()
indirect_predictor= 0             #Indirect predictor included? 0 for no; 1 for 1 level upper; 2 for all possible ones

result_outfiles   = c('result_matrix.csv','outer_q2.csv','beta_matrix.csv')

# Internal parameters
cell_col          = 11:20
total_cell_col    = 9:18
num_test          = 1 #number of test sets for outer loop

# Functions:
# given a site_id, return a list with list$x are the predictors, and list$y are true responses
data_prepare <- function(site_id,msdata,totdata,network,prot_interact = prot_interact, gene_interact = gene_interact,
                         cell_col=cell_col,total_cell_col=total_cell_col,mutation_matrix = mutation_matrix){
  data_model <- list(x=c(),y=c())
  protein_name <- msdata[site_id,"merge_id"]
  gene_name <- msdata[site_id,'gene_symbol']
  site_list <- strsplit(msdata[site_id,'Site.Position'],split = ';')[[1]]
  all_predictors_genes <- names(network[network[,gene_name] > 0,gene_name])
  all_predictors_data  <- msdata[msdata$gene_symbol %in% all_predictors_genes,]
  if(nrow(all_predictors_data) < 1){
    return(NULL)
  }else{
    relate_evidence <- prot_interact[prot_interact$protB == protein_name,]
    site_evidence   <- relate_evidence[relate_evidence$site %in% site_list,'protA']
    prot_evidence   <- relate_evidence[!(prot_interact$site[prot_interact$protB == protein_name] %in% site_list),'protA']
    gene_evidence   <- gene_interact[gene_interact$geneA == gene_name | gene_interact$geneB == protein_name,]
    if(nrow(gene_evidence) > 0){
      gene_evidence  <- unlist(gene_evidence[gene_evidence != gene_name])
    }
    penalty_vector  <- rep(3,nrow(all_predictors_data))
    penalty_vector[which(all_predictors_data$gene_symbol %in% gene_evidence)] <- 2
    penalty_vector[which(all_predictors_data$merge_id %in% prot_evidence)]     <- 1
    penalty_vector[which(all_predictors_data$merge_id %in% site_evidence)]     <- 0.001
    tot_protein <- as.numeric(totdata[totdata$acc_id == strsplit(site_id,split = '_')[[1]][1],total_cell_col])
    data_model$x <- t(rbind(all_predictors_data[,cell_col],'tot_protein' = tot_protein, 'mutation' = as.numeric(mutation_matrix[gene_name,])))
    data_model$y <- as.numeric(msdata[site_id,cell_col])
    data_model$penalty <- c(penalty_vector,1,1)
  }
  return (data_model)
}

# Nested cross validation: alpha and lambda are vectors for grid search, if lambda is set to NULL, it would train alpha only and let glmnet to decide lambda
# inner loop would use a leave-1-out CV, outer loop would have ntest for test set.
nestcv <- function(model_data,alphas,lambdas,ntest = num_test){
  cell_line_list <- rownames(model_data$x)
  penalty <- model_data$penalty
  result_matrix <- matrix(nrow = 0,ncol=5)
  colnames(result_matrix) <- c('cell_out','alpha','inner_q2','predict_value','true_value')
  outer_pred_list <- rep(0,nrow(model_data$x))
  beta_matrix <- matrix(0,nrow = 0,ncol = ncol(model_data$x)+1) #include all coefficient and intercept
  ind <- c(1:10,1:10)
  beta0 <- c()
  # outer loop
  for (i in 1:10){
    cell_out <- cell_line_list[i]
    test <- ind[i:(i+ntest-1)]
    outer_test_x <- model_data$x[test,]
    if(ntest == 1){outer_test_x <- t(outer_test_x)}
    outer_train_x <- model_data$x[-test,]
    outer_test_y <- model_data$y[test]
    outer_train_y <- model_data$y[-test]
    best_alpha <- 0
    best_inner_q2 <- -Inf
    # inner loop
    if(is.null(lambdas)){
      # lambda not given, train alpha only
      for (alpha in alphas){
        best_lambda <- NULL
        inner_predicts <- lasso_cv(x = outer_train_x,y = outer_train_y,alpha = alpha, lambda = NULL,ntest = 1,penalty = penalty)
        q2 <- 1 - sum((inner_predicts-outer_train_y)^2)/sum((outer_train_y-mean(outer_train_y))^2)
        if(q2 > best_inner_q2){
          best_alpha <- alpha
          best_inner_q2 <- q2
        }
      }
    }else{
      #grid search on both lambda and alpha
      for (lambda in lambdas){
        for (alpha in alphas){
          inner_predicts <- lasso_cv(x = outer_train_x,y = outer_train_y,alpha = alpha, lambda = lambda,penalty = penalty)
          q2 <- 1 - sum((inner_predicts-outer_train_y)^2)/sum((outer_train_y-mean(outer_train_y))^2)
          if(q2 > best_inner_q2){
            best_alpha <- alpha
            best_lambda <- lambda
            best_inner_q2 <- q2
          }
        }
      }
    }
    #end on inner loop
    # use best alpha and lambda on test data
    outer_test_model <- glmnet(outer_train_x,outer_train_y,alpha = best_alpha,lambda = best_lambda,penalty.factor = penalty)
    return_lambda <- outer_test_model$lambda[which.max(outer_test_model$dev.ratio)]
    if(is.nan(return_lambda) | is.infinite(return_lambda)){
      outer_test_model <- glmnet(outer_train_x,outer_train_y,alpha = best_alpha,lambda = 0.001,penalty.factor = penalty)
    }
    outer_test_pred <- predict(outer_test_model,newx = outer_test_x,s = return_lambda)
    if(is.null(best_lambda)){best_lambda <- return_lambda}
    result_matrix <- rbind(result_matrix,c(cell_out,best_alpha,best_inner_q2,outer_test_pred,outer_test_y))
    beta0 <- outer_test_model$a0[which.max(outer_test_model$dev.ratio)]
    beta_matrix <- rbind(beta_matrix,c('intercept' = as.numeric(beta0),outer_test_model$beta[,which.max(outer_test_model$dev.ratio)]))
  }
  best_outer_q2 <- 1 - sum((as.numeric(result_matrix[,'predict_value'])-model_data$y)^2)/sum((model_data$y-mean(model_data$y))^2)
  return(list('matrix' = result_matrix,'best_outer_q2' = best_outer_q2,'beta_matrix' = beta_matrix))
}

# leave-1-out cross validation,return predicted values
lasso_cv <- function(x,y,alpha,lambda = NULL,ntest = 1,penalty = NULL){
  pred <- rep(0,nrow(x))
  for(i in 1:(nrow(x) +1 -ntest)){
    test <- i:(i+ntest-1)
    fit_cv <- glmnet(x[-test,],y[-test],lambda = lambda,alpha = alpha,penalty.factor = penalty)
    best_lambda <- fit_cv$lambda[which.max(fit_cv$dev.ratio)]
    if(is.nan(best_lambda) | is.infinite(best_lambda)){
      fit_cv <- glmnet(x[-test,],y[-test],lambda = 0.001,alpha = alpha,penalty.factor = penalty)
    }
    pred[i] <- predict(fit_cv,t(x[test,]),s = best_lambda)
  }
  return(pred)
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

# Main body
msdata_raw      <- read.csv(paste(script_dir,msdata_filename,sep = '/'),header = T,as.is = T)
totdata         <- read.csv(paste(script_dir,tot_prt_filename,sep = '/'),header = T,as.is = T)
network         <- as.matrix(read.csv(paste(network_dir,network_filename,sep = '/'),row.names = 1))
mutation_matrix <- read.csv(paste(script_dir,mutation_filename,sep = '/'),header = T,as.is = T,row.names = 1)
prot_interact   <- read.csv(paste(network_dir,prot_table_file,sep = '/'),as.is = T)
gene_interact   <- read.csv(paste(network_dir,gene_table_file,sep = '/'),as.is = T)
msdata          <- msdata_raw[msdata_raw$gene_symbol %in% rownames(network) & msdata_raw$merge_id %in% totdata$merge_id,]
unique_id       <- paste(msdata$acc_id,msdata$Site.Position,sep = '_')
gene_list       <- msdata$gene_symbol
rownames(msdata) <- unique_id


result_matrix <- matrix(0,ncol = 7,nrow = 0)
outerq_matrix <- matrix(0,ncol = 3,nrow = 0)
beta_matrix   <- matrix(0,ncol=length(cell_col)+3,nrow = 0)

colnames(result_matrix) <- c('gene_symbol','siteid','cell_out','best_alpha','best_inner_q2','best_pred','true_value')
colnames(outerq_matrix) <- c('site_id','gene','outer_q2')
colnames(beta_matrix)   <- c('site_id','gene','predictor',colnames(msdata[cell_col]))

# iterate of all phos sites related to MAPK pathway
for (i in 1:length(unique_id)){
  site_id     <- unique_id[i]
  gene_symbol <- gene_list[i]
  model_data  <- data_prepare(site_id,msdata,totdata = totdata,network = network,prot_interact = prot_interact,
                              gene_interact = gene_interact,cell_col=cell_col,total_cell_col = total_cell_col,mutation_matrix = mutation_matrix)
  name_y <- paste(gene_symbol,site_id,sep = '_')
  if (!is.null(model_data)){
    model_result  <- nestcv(model_data,alphas = alphas,lambdas = lambdas,ntest = num_test)
    new_result    <- cbind('gene' = rep(gene_symbol,length(cell_col)),'site_id' = rep(site_id,length(cell_col)),model_result$matrix)
    result_matrix <- rbind(result_matrix,new_result)
    outerq_matrix <- rbind(outerq_matrix,cbind(site_id,gene_symbol,model_result$best_outer_q2))
    beta_matrix   <- rbind(beta_matrix,cbind(site_id,gene_symbol,colnames(model_result$beta_matrix),t(model_result$beta_matrix)))
  }
}

write.csv(result_matrix,paste(results_dir,result_outfiles[1],sep = '/'),row.names = F)
write.csv(outerq_matrix,paste(results_dir,result_outfiles[2],sep = '/'),row.names = F)
write.csv(beta_matrix,paste(results_dir,result_outfiles[3],sep = '/'),row.names = F)