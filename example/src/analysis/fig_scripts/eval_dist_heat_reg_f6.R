dist_range= 1:5
configs <- paste(c('parallel_config_en_large_dist_k'),dist_range,sep = '')
heat_k <- c('1','01','001','0008','0006','0004','0002','0001','00001','000001')

configs_heat <- paste('parallel_config_en_large_heat_ud_k',heat_k,sep = '')
output_path <- '~/Documents/workspace/phospho_network/example/script_files/output/'
en_all_path <- 'enL_all_results'
random_table <- read.csv('~/Documents/workspace/phospho_network/example/script_files/result_integrate/random_integrate.csv',as.is = T)

result_table <- function(result_path){
  all_files    <- list.files(result_path)
  result_files <- grep('result',all_files,value = T)
  result       <- NULL
  for(result_file in result_files){
    if(is.null(result)){
      result <- read.csv(paste(result_path,result_file,sep = '/'),as.is = T)
    }else{
      result <- rbind(result,read.csv(paste(result_path,result_file,sep = '/'),as.is = T))
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
best_regression <- function(configs,output_path){
  best_result <- NULL
  for (i in 1:length(configs)){
    result <- result_table(paste(output_path,configs[i],sep = ''))
    result <- result[result$best_outer_q2>0,]
    if(is.null(best_result) & nrow(result)>0){
      best_result <- cbind(result,'best_config'=configs[i])
    }else{
      all_targets <- unique(result$gene_site)
      for(target in all_targets){
        if(!target %in% best_result$gene_site){
          best_result <- rbind(best_result,cbind(result[result$gene_site==target,],'best_config'=configs[i]))
        }else if(unique(best_result$best_outer_q2[best_result$gene_site==target]) < unique(result$best_outer_q2[result$gene_site==target])){
          best_result <- rbind(best_result[best_result$gene_site != target,],cbind(result[result$gene_site==target,],'best_config'=configs[i]))
        }
      }
    }
  }
  return(best_result)
}
  
dist_analysis <- function(configs,output_path,target_limit){
  dist_result <- NULL
  for (i in 1:length(configs)){
    best_result <- best_regression(configs[1:i],output_path)
    best_result <- best_result[best_result$gene_site %in% target_limit,]
    dist_result <- rbind(dist_result,c('best_q2' = cal_q2(best_result$predict_value,best_result$true_value),'median_q2' = median(best_result$best_outer_q2),'best_rmse'=cal_rmse(best_result$predict_value,best_result$true_value),'dist_scope' = configs[i],'mpred_given' = median(best_result$predictors_given)))
  }
  dist_result <- data.frame(dist_result,stringsAsFactors = F)
  dist_result$mpred_given <- as.numeric(dist_result$mpred_given)
  dist_result$best_q2 <- as.numeric(dist_result$best_q2)
  dist_result$median_q2 <- as.numeric(dist_result$median_q2)
  dist_result$best_rmse <- as.numeric(dist_result$best_rmse)
  return(dist_result)
}

target_limit <- Reduce(intersect,list(best_regression(configs[1:2],output_path)[,1],best_regression(configs_heat[2],output_path)[,1]))
dist_result <- dist_analysis(configs[-1],output_path,target_limit)
heat_result <- dist_analysis(configs_heat,output_path,target_limit)

par(mfrow=c(1,2))
plot(log(dist_result$mpred_given,base = 10),dist_result$median_q2,
     xlim = log(c(min(dist_result$mpred_given,heat_result$mpred_given),max(dist_result$mpred_given,heat_result$mpred_given)),base = 10),
     ylim = c(min(dist_result$median_q2,heat_result$median_q2,random_table[,1]),max(dist_result$median_q2,heat_result$median_q2)),
     xlab = expression('Log'[10]*' median number of predictors'),ylab = 'Median Q-squared',main = 'Q-squared vs. Predictor number')
lines(log(dist_result$mpred_given,base = 10),dist_result$median_q2,col='blue')
points(log(heat_result$mpred_given,base = 10),heat_result$median_q2)
lines(log(heat_result$mpred_given,base = 10),heat_result$median_q2,col='red')
points(log(random_table[,5],base = 10),random_table[,1])
lines(log(random_table[,5],base = 10),random_table[,1])
legend('topleft',legend = c('heat diffussion','distance','Random Select'),col = c('red','blue','black'),lty = 1,bty = 'n',cex = 0.8)
abline(v = log(heat_result$mpred_given[8],10),lty=2)

plot(log(dist_result$mpred_given,base = 10),dist_result$best_rmse,
     xlim = log(c(min(dist_result$mpred_given,heat_result$mpred_given),max(dist_result$mpred_given,heat_result$mpred_given)),base = 10),
     ylim = c(min(dist_result$best_rmse,heat_result$best_rmse),max(dist_result$best_rmse,heat_result$best_rmse,random_table[,3])),
     xlab = expression('Log'[10]*' median number of predictors'),ylab = 'Median RMSE',main = 'RMSE vs. Predictor number')
lines(log(dist_result$mpred_given,base = 10),dist_result$best_rmse,col='blue')
points(log(heat_result$mpred_given,base = 10),heat_result$best_rmse)
lines(log(heat_result$mpred_given,base = 10),heat_result$best_rmse,col='red')
points(log(random_table[,5],base = 10),random_table[,3])
lines(log(random_table[,5],base = 10),random_table[,3])
legend(3.11,0.15847,legend = c('heat diffussion','distance','Random Select'),col = c('red','blue','black'),lty = 1,bty = 'n',cex = 0.8)
abline(v = log(heat_result$mpred_given[8],10),lty=2)
