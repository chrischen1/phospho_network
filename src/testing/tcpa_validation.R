# TCPA validation
# data parser and analysis

library(glmnet)
library(methods)
library(igraph)

tcpa_data_file        = '~/Documents/workspace/phospho_network/processed_data/tcpa/tcpa_data_processed.csv'
interaction_site_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/interaction_table_site_all.csv'
interaction_gene_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/network_all_gene.csv'
results_dir           = '~/Documents/workspace/phospho_network/script_files/analysis_rppa'
result_outfiles       = c('result_matrix.csv','outer_q2.csv','beta_matrix.csv')

data_col <- -(1:3)
penalty_site      = 0.01
penalty_prot      = 1
penalty_gene      = 2


expand_matrix <- function(m,ncol,sep = '; '){
  m2 <- m[0,]
  for(i in 1:nrow(m)){
    var_split <- strsplit(m[i,ncol],split = sep)[[1]]
    if(length(var_split) <= 1){
      m2 <- rbind(m2,m[i,])
    }else if(length(var_split) > 1){
      new_table <- matrix(rep(as.character(m[i,]),length(var_split)),nrow = length(var_split),byrow = T)
      new_table[,ncol] <- var_split
      colnames(new_table) <- colnames(m2)
      m2 <- rbind(m2,new_table)
    }
  }
  return(m2)
}

data_prepare <- function(site_id,msdata,totdata,network,prot_interact = prot_interact, gene_interact = NULL,
                         data_col=data_col,total_data_col=data_col,mutation_matrix = NULL,penalty_site= 0.01,penalty_prot= 1,penalty_gene= 2,penalty_noevd = 3){
  data_model <- list(x=c(),y=c())
  gene_name <- strsplit(site_id,split = '_')[[1]][1]
  site_list <- strsplit(msdata[site_id,'site'],split = ';')[[1]]
  all_predictors_genes <- rownames(network)[network[,gene_name] > 0]
  all_predictors_data  <- msdata[msdata$gene_symbol %in% all_predictors_genes,]
  if(nrow(all_predictors_data) < 1){
    return(NULL)
  }else{
    penalty_vector  <- rep(penalty_noevd,nrow(all_predictors_data))

    if(!is.null(prot_interact)){
      relate_evidence <- prot_interact[prot_interact$geneB == gene_name,]
      prot_evidence   <- relate_evidence[!(prot_interact$site[prot_interact$geneB == gene_name] %in% site_list),'geneA']
      site_evidence   <- relate_evidence[relate_evidence$site %in% site_list,'geneA']
      penalty_vector[which(all_predictors_data$gene_symbol %in% prot_evidence)]     <- penalty_prot
      penalty_vector[which(all_predictors_data$gene_symbol %in% site_evidence)]     <- penalty_site
    }
    if (!is.null(gene_interact)){
      gene_evidence   <- gene_interact[gene_interact$geneB == gene_name,]
      if(nrow(gene_evidence) > 0){
        gene_evidence  <- unlist(gene_evidence[gene_evidence != gene_name])
        penalty_vector[which(all_predictors_data$gene_symbol %in% gene_evidence)]  <- penalty_gene
      }
    }
    x <- t(all_predictors_data[,data_col])
    if(gene_name %in% totdata$gene_symbol){
      tot_protein <- as.numeric(totdata[totdata[,'gene_symbol'] == gene_name,total_data_col])
      x <- cbind(x,'tot_protein' = tot_protein)
      penalty_vector <- c(penalty_vector,1)
    }
    if(!is.null(mutation_matrix)){
      x <- cbind(x,'mutation' = as.numeric(mutation_matrix[gene_name,]))
      penalty_vector <- c(penalty_vector,1)
    }
    data_model$x <- x
    data_model$y <- as.numeric(msdata[site_id,data_col])
    data_model$penalty <- penalty_vector
  }
  return (data_model)
}

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

# main body
tcpa_data_raw <- read.csv(tcpa_data_file,as.is = T,header = T)
interaction_site_raw <- read.csv(interaction_site_file,as.is = T)
interaction_gene_raw <- read.csv(interaction_gene_file,as.is = T)

tcpa_data     <- expand_matrix(tcpa_data_raw,ncol = which(colnames(tcpa_data_raw) == 'gene_symbol'))
tcpa_data_flt <- tcpa_data[tcpa_data[,'gene_symbol'] %in% unique(unlist(interaction_gene_raw)),]
tcpa_data_tot <- tcpa_data_flt[tcpa_data_flt[,'site'] == '',]
tcpa_data_pho <- tcpa_data_flt[tcpa_data_flt[,'site'] != '',]
rownames(tcpa_data_pho) <- paste(tcpa_data_pho[,'gene_symbol'],tcpa_data_pho[,'site'],sep = '_')
tcpa_genes    <- unique(tcpa_data_flt[,'gene_symbol'])

prot_interact    <- interaction_site_raw[interaction_site_raw$geneA %in% tcpa_genes & interaction_site_raw$geneB %in% tcpa_genes & interaction_site_raw[,'geneA'] != interaction_site_raw[,'geneB'],]
interaction_gene <- interaction_gene_raw[interaction_gene_raw$geneA %in% tcpa_genes & interaction_gene_raw$geneB %in% tcpa_genes & interaction_gene_raw$geneA != interaction_gene_raw$geneB,]
g <- graph_from_edgelist(as.matrix(interaction_gene))
network <- as.matrix(as_adjacency_matrix(g))

# iterate of all phos sites on the dataset

for (i in 1:nrow(tcpa_data_pho)){
  site_id     <- rownames(tcpa_data_pho)[i]
  gene_symbol <- tcpa_data_pho[i,'gene_symbol']
  print(site_id)
  model_data  <- data_prepare(site_id,tcpa_data_pho,totdata = tcpa_data_tot,network = network,prot_interact = prot_interact,
                              data_col=data_col, total_data_col=data_col,penalty_site = penalty_site,penalty_prot = penalty_prot,penalty_gene = penalty_gene)
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