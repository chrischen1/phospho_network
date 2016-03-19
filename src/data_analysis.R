# This Rscript would take a processed MS data file(data_processing.R) as input, and do the following steps:
# 1.Import the MS data and network file(as an adjacency matrix) from csv files
# 2.Filter the data by keeping only proteins in the network 
# 3.Nodes 1 level further upstream would also be treated as predictors, but a penalty factor would be imposed
# 4.Construct the LASSO model: PeptideA_site1 ~ PeptideA_site1 + PeptideA_site2...
# 5.Evaluate the model with r-square and q-square

library(glmnet)

# Set parameters
script_dir        = '~/Documents/workspace/phospho_network/script_files'
results_dir       = '~/Documents/workspace/phospho_network/script_files/analysis_results'
network_dir       = '~/Documents/workspace/phospho_network/script_files/interactions_mapk'
msdata_filename   = 'msdata_processed.csv'
tot_prt_filename  = 'total_protein_processed.csv'
network_filename  = 'mapk_network.csv'
evidence_filename = ''
mutation_filename = 'mutation_matrix.csv'
alphas            = seq(0,1,0.05)
lambdas           = c()
sig_level         = 0.1
indirect_predictor= 1   #Indirect predictor included? 0 for no; 1 for 1 level upper; 2 for all possible ones
penalty_sequence  = seq(1,10,0.5)

# Internal parameters
cell_col          = 11:20
total_cell_col    = 9:18
num_test          = 1 #number of test sets for outer loop

# Functions:
# given a site_id, return a list with list$x are the predictors, and list$y are true responses
data_prepare <- function(site_id,msdata,totdata,network,cell_col=cell_col,total_cell_col=total_cell_col,mutation_matrix = mutation_matrix){
  data_model <- list(x=c(),y=c())
  protein_name <- msdata[site_id,"merge_id"]
  gene_name <- msdata[site_id,'gene_symbol']
  predictors_names <- rownames(network)[network[,colnames(network) == protein_name] > 0]
  data_predictor <- msdata[msdata$merge_id %in% predictors_names,]
  if (nrow(data_predictor) < 2){
    return(NULL)
  }else{
    tot_protein <- as.numeric(totdata[totdata$acc_id == strsplit(site_id,split = '_')[[1]][1],total_cell_col])
    data_model$x <- t(rbind(data_predictor[,cell_col],'tot_protein' = tot_protein, 'mutation' = as.numeric(mutation_matrix[gene_name,])))
    data_model$y <- as.numeric(msdata[site_id,cell_col])
    data_model$penalty <- c(network[data_predictor$merge_id,protein_name],'tot_protein'=1,'mutation'=1)
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
msdata <- read.csv(paste(script_dir,msdata_filename,sep = '/'),header = T,as.is = T)
totdata <- read.csv(paste(script_dir,tot_prt_filename,sep = '/'),header = T,as.is = T)
network <- read.csv(paste(network_dir,network_filename,sep = '/'))[,-1]
tot_merge_id <- gsub('(\\w+)-\\d+','\\1',totdata$acc_id)
totdata[,1] <- tot_merge_id
colnames(totdata)[1] <- 'merge_id'
rownames(network) <- colnames(network)
network <- as.matrix(network)

mutation_matrix_raw <- read.csv(paste(script_dir,mutation_filename,sep = '/'),header = T,as.is = T)
mutation_matrix <- mutation_matrix_raw[,-1]
rownames(mutation_matrix) <- mutation_matrix_raw[,1]

msdata_filtered <- msdata[msdata$merge_id %in% rownames(network) & msdata$merge_id %in% totdata$merge_id,]
unique_id <- paste(msdata_filtered$acc_id,msdata_filtered$Site.Position,sep = '_')
rownames(msdata_filtered) <- unique_id
gene_list <- msdata_filtered[unique_id,]$gene_symbol

node_matrix <- matrix(0,nrow = 0,ncol = 4)
colnames(node_matrix) <- c('node','label','outer_q2','p_val')
beta_matrix <- matrix(0,nrow = 0,ncol = 6)
colnames(beta_matrix) <- c('nodeA','nodeB','outer_q2','p_val','cv_mbeta','overal_beta')
edge_matrix <- matrix(0,nrow = 0,ncol = 6)
colnames(edge_matrix) <- c('nodeA','nodeB','abs_mbeta','sign_mbeta','outer_q2','p_val')
result_matrix <- matrix(0,ncol = 9,nrow = 0)
colnames(result_matrix) <- c('gene_symbol','siteid','outer_q2','p_val','cell_out','best_alpha','best_inner_q2','best_pred','true_value')
error_matrix <- matrix(0,ncol=length(cell_col)+1,nrow = 0)
colnames(error_matrix) <- c('gene',colnames(msdata_filtered[cell_col]))

# iterate of all phos sites related to MAPK pathway
penalty_sequence <- c(penalty_sequence,1000) # A very large penalty would be the same as prediction without indirect nodes
par(mfrow=c(1,2))
matrix_1up <- matrix(0,nrow = 0,ncol = length(penalty_sequence))
colnames(matrix_1up) = penalty_sequence
matrix_aup <- matrix_1up
names_y <- c()
q2_vanila <- c()

for (i in 1:length(unique_id)){
  site_id <- unique_id[i]
  gene_symbol <- gene_list[i]
  model_data <- data_prepare(site_id,msdata_filtered,totdata = totdata,network = network,cell_col=cell_col,total_cell_col = total_cell_col,mutation_matrix = mutation_matrix)
  name_y <- paste(gene_symbol,site_id,sep = '_')
  sample_num <- length(model_data$y)
  if (!is.null(model_data)){
    names_y <- c(names_y,name_y)
    model_vanila <- nestcv(model_data,alphas = alphas,lambdas = lambdas,ntest = num_test)
    q2_vanila<- c(q2_vanila,model_vanila$best_outer_q2)
    q2_1up <- c();q2_aup <- c()
    for (j in penalty_sequence){
      model_data_1up <- data_prepare(site_id,msdata_filtered,totdata = totdata,network = up_network(network,all_level = F,penalty = j),
                                     cell_col=cell_col,total_cell_col = total_cell_col,mutation_matrix = mutation_matrix)
      model_1up <- nestcv(model_data_1up,alphas = alphas,lambdas = lambdas,ntest = num_test)
      q2_1up <- c(q2_1up,model_1up$best_outer_q2)
      model_data_aup <- data_prepare(site_id,msdata_filtered,totdata = totdata,network = up_network(network,all_level = T,penalty = j),
                                     cell_col=cell_col,total_cell_col = total_cell_col,mutation_matrix = mutation_matrix)
      model_aup <- nestcv(model_data_aup,alphas = alphas,lambdas = lambdas,ntest = num_test)
      q2_aup <- c(q2_aup,model_aup$best_outer_q2)
    }
    matrix_1up <- rbind(matrix_1up,q2_1up)
    matrix_aup <- rbind(matrix_aup,q2_aup)
  }
}
matrix_1up <- cbind(matrix_1up,q2_vanila)
matrix_aup <- cbind(matrix_aup,q2_vanila)
write.csv(matrix_1up,file = paste(results_dir,'matrix_1up.csv',sep = '/'),row.names = names_y)
write.csv(matrix_aup,file = paste(results_dir,'matrix_aup.csv',sep = '/'),row.names = names_y)



