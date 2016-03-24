# TCPA validation
# data parser and analysis

library(glmnet)
library(methods)
library(igraph)

tcpa_data_file        = '~/Documents/workspace/phospho_network/processed_data/tcpa/tcpa_data_processed.csv'
interaction_site_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/interaction_table_site_all.csv'
interaction_gene_file = '~/Documents/workspace/phospho_network/script_files/interaction_data_processed/network_all_gene.csv'
target_gene_file      = '~/Documents/workspace/phospho_network/script_files/analysis_rppa/response_list.txt'

results_dir           = '~/Documents/workspace/phospho_network/script_files/analysis_rppa'
result_outfiles       = c('result_matrix.csv','beta_matrix.csv')

data_col          = 4:413
penalty_site      = 0.01
penalty_prot      = 1
penalty_gene      = 2
alphas            = seq(0,1,0.05)
lambdas           = c()
indirect_predictor= 0
num_test          = 1 #number of test sets for outer loop
outerfold         = 3

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

data_prepare <- function(site_id,msdata,totdata,network,prot_interact = prot_interact, gene_interact = NULL,
                         data_col=data_col,total_data_col=data_col,mutation_matrix = NULL,penalty_site= 0.01,penalty_prot= 1,penalty_gene= 2,penalty_noevd = 3){
  data_model <- list(x=c(),y=c())
  gene_name <- strsplit(site_id,split = '_')[[1]][1]
  site_list <- strsplit(msdata[site_id,'site'],split = ';')[[1]]
  all_predictors_genes <- rownames(network)[network[,gene_name] > 0]
  all_predictors_data  <- msdata[msdata$gene_symbol %in% all_predictors_genes,]
  if('collision' %in% colnames(msdata)){
    if(msdata[site_id,'collision'] != ''){
      collision_names <- strsplit(msdata[site_id,'collision'],'; ')[[1]]
      collision_ids   <- paste(collision_names,msdata[site_id,'site'],sep = '_')
      all_predictors_data <- all_predictors_data[-(which(rownames(all_predictors_data) %in% collision_ids)),]
    }
  }
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
    x2 <- apply(x,2,as.numeric)
    if(ncol(x2) < 2){return (NULL)}
    rownames(x2) <- rownames(x)
    data_model$x <- x2
    data_model$y <- as.numeric(msdata[site_id,data_col])
    data_model$penalty <- penalty_vector
    data_model$site_id <- site_id
  }
  return (data_model)
}

nestcv <- function(model_data,alphas,lambdas = NULL,outerfold = 3,innerfold = 10){
  penalty         <- model_data$penalty
  test_list       <- list()
  sample_list     <- 1:nrow(model_data$x)
  sample_list_tmp <- sample_list
  for (i in 1:outerfold){
    test_ind <- sample(sample_list_tmp,size = floor(length(sample_list)/outerfold),replace = F)
    test_list[[i]] <- test_ind
    sample_list_tmp <- sample_list_tmp[-test_ind]
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
      cv_model <- cv.glmnet(model_data$x[-test_list[[i]],],model_data$y[-test_list[[i]]],nfolds = innerfold,penalty.factor = model_data$penalty,alpha = alpha)
      return_err <- min(cv_model$cvm)
      if(return_err < best_err){
        best_alpha <- alpha
        best_model <- cv_model$glmnet.fit
      }
    }
    #end on inner loop
    # use best alpha on test data
    outer_test_pred <- predict(best_model,newx = model_data$x[outer_test_set,],s = cv_model$lambda.min)
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
tar_genes  <- gsub('(\\w+)_.+','\\1',unlist(read.table(target_gene_file,as.is = T)))

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
result_matrix <- matrix(0,nrow = 0,ncol = 6)
beta_matrix   <- matrix(0,nrow = 0,ncol = 2+outerfold)
colnames(beta_matrix)   <- c('gene_site','predictor',paste('test',1:outerfold,sep = '_'))
colnames(result_matrix) <- c('gene_site', 'test_set', 'alpha', 'predict_value', 'true_value', 'best_outer_q2')

target_list <- rownames(tcpa_data_pho)[tcpa_data_pho$gene_symbol %in% tar_genes]

for (i in 1:length(target_list)){
  site_id     <- target_list[i]
  gene_symbol <- strsplit(site_id,split = '_')[[1]][1]
  model_data  <- data_prepare(site_id,tcpa_data_pho,totdata = tcpa_data_tot,network = network,prot_interact = prot_interact,
                              data_col=data_col, total_data_col=data_col,penalty_site = penalty_site,penalty_prot = penalty_prot,penalty_gene = penalty_gene)
  name_y <- paste(gene_symbol,site_id,sep = '_')
  if (!is.null(model_data)){
    model_result  <- nestcv(model_data,alphas = alphas,outerfold = outerfold)
    result_matrix <- rbind(result_matrix,model_result$matrix)
    beta_matrix   <- rbind(beta_matrix,model_result$beta_matrix)
  }
}

write.csv(result_matrix,paste(results_dir,result_outfiles[1],sep = '/'),row.names = F)
write.csv(beta_matrix,paste(results_dir,result_outfiles[2],sep = '/'),row.names = F)
