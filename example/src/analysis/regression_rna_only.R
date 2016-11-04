mdata_filename = '~/Documents/workspace/phospho_network/example/script_files/rppa_processed.csv'
rna_filename   = '~/Documents/workspace/phospho_network/example/script_files/rna_processed.csv'
out_filename   = '~/Documents/workspace/phospho_network/example/script_files/output/pred_rna.csv'
nfold          = 10




sample_prepare <- function(nsample,fold){
  test_list       <- list()
  train_list      <- list()
  sample_list     <- 1:nsample
  sample_list_tmp <- sample(sample_list,size = floor(length(sample_list)/fold)*fold)
  test_size <- floor(nsample/fold)
  for (i in 1:fold){
    test_list[[i]] <- sample_list_tmp[1:test_size]
    sample_list_tmp <- sample_list_tmp[-(1:test_size)]
    train_list[[i]] <- sample_list[-test_list[[i]]]
  }
  return(list('train' = train_list,'test' = test_list))
}

cv_lm <- function(target_val,pred_val,nfold){
  pred_table <- NULL
  sample_list <- sample_prepare(length(target_val),nfold)
  for(i in 1:nfold){
    train_x <- pred_val[sample_list$train[[i]],]
    train_y <- as.numeric(target_val[sample_list$train[[i]]])
    test_x <- data.frame(pred_val[sample_list$test[[i]],])
    test_y <- as.numeric(target_val[sample_list$test[[i]]])
    pred_model <- lm(train_y~train_x)
    pred_y <- as.numeric(pred_model$coefficients %*%  t(cbind(1,test_x)))
    pred_table <- rbind(pred_table,cbind('true' = test_y,'pred' = pred_y))
  }
  return(pred_table)
}

data_norm <- function(x){
  if(max(x) != min(x)){return ((x-min(x))/(max(x)-min(x)))}
  else{return(rep(1,length(x)))}
}

rm_na <- function(x){
  x[is.na(x)] <- mean(x[!is.na(x)])
  return(x)
}

mdata     <- read.csv(mdata_filename,row.names = 1)
rna_data  <- read.csv(rna_filename,as.is = T,row.names = 1)
target_list <- rownames(mdata)[grep('_\\w+',rownames(mdata))]
target_genes <- strsplit(unlist(lapply(strsplit(target_list,'_'),function(x)x[1])),';')

mdata <-apply(mdata,2,data_norm)
rna_data <- apply(rna_data,2,function(x)data_norm(rm_na(x)) )

pred_table_all <- NULL
for(i in 1:length(target_genes)){
  print(i)
  pred_gene <- target_genes[[i]]
  target <- target_list[i]
  target_val <- mdata[target,]
  if(sum(rownames(rna_data) %in% pred_gene) > 1){
    pred_val <- t(rna_data[rownames(rna_data) %in% pred_gene,])
    pred_table_all <- rbind(pred_table_all,cbind('target'=target,cv_lm(target_val,pred_val,nfold)))
  }else if(sum(rownames(rna_data) %in% pred_gene) == 1){
    pred_val <- t(t(rna_data[rownames(rna_data) %in% pred_gene,]))
    pred_table_all <- rbind(pred_table_all,cbind('target'=target,cv_lm(target_val,pred_val,nfold)))
  }
}

write.csv(pred_table_all,out_filename,row.names = F)









