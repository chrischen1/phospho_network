random_path = '~/Documents/workspace/phospho_network/example/script_files/output/random'
random_configs = c('enL_random_results_k100','enL_random_results_k300','enL_random_results_k900','enL_random_results_k2700','enL_random_results_k8100')
output_file    = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/random_integrate.csv'

get_result <- function(config_path,keyword){
  all_files    <- list.files(config_path)
  result_files <- grep(keyword,all_files,value = T)
  result       <- NULL
  for(result_file in result_files){
    if(is.null(result)){
      result <- read.csv(paste(config_path,result_file,sep = '/'),as.is = T)
    }else{
      result <- rbind(result,read.csv(paste(config_path,result_file,sep = '/'),as.is = T))
    }
  }
  return(result)
}
cal_q2 <- function(x,y){
  return(round(1-sum((x-y)^2)/sum((y-mean(y))^2),4))
}
cal_rmse <- function(x,y){
  return(sqrt(sum((x-y)^2)/length(x)))
}
random_table <- NULL
for (config in random_configs){
  result_config <- get_result(paste(random_path,config,sep = '/'),keyword = 'result')
  for(target in letters[1:14]){
    target_q2 <- c()
    target_rmse <- c()
    for(i in 1:20){
      target_q2 <- c(target_q2,unique(result_config[result_config[,1] == paste(target,i,sep = '_'),'best_outer_q2']))
      target_rmse <- c(target_rmse,cal_rmse(result_config[result_config[,1] == paste(target,i,sep = '_'),4],
                                            result_config[result_config[,1] == paste(target,i,sep = '_'),5]))
    }
    target_q2_median <- median(target_q2)
    target_q2_best <- max(target_q2)
    target_rmse_median <- median(target_rmse)
    target_rmse_best <- min(target_rmse)
    pred_num <- unique(result_config[,7])
    random_table <- rbind(random_table,cbind('target'=target,target_q2_median,target_q2_best,target_rmse_median,target_rmse_best,pred_num))
  }
}
random_table_all <- NULL
for(i in unique(random_table[,6])){
  config_table <- random_table[random_table[,6]==i,]
  random_table_all <- rbind(random_table_all,c(apply(config_table[,2:5],2,function(x)median(as.numeric(x))),'pred_num'=i))
}

write.csv(random_table_all,output_file,row.names = F)
