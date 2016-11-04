dist_range= 1:4
configs <- paste(c('parallel_config_en_large_dist_k'),dist_range,sep = '')
heat_k <- c('1','01','001','0008','0006','0004','0002','0001','00001','000001')

configs_heat <- paste('parallel_config_en_large_heat_ud_k',heat_k,sep = '')
output_path <- '~/Documents/workspace/phospho_network/example/script_files/output/'
en_all_path <- 'enL_all_results'
  
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
dist_analysis <- function(configs,output_path){
  best_q2 <- c()
  best_k <- c()
  best_pred_num <-c()
  best_rmse <- c()
  for (i in 1:length(configs)){
    result <- result_table(paste(output_path,configs[i],sep = ''))
    result <- result[result$best_outer_q2>0,]
    if(nrow(result)>0){
      all_target <- unique(result$gene_site)
      for(target in all_target){
        new_q2 <- unique(result$best_outer_q2[result$gene_site == target])
        new_k <- configs[i]
        new_pnum <- unique(result$predictors_given[result$gene_site == target])
        new_rmse <- sqrt(sum((result$predict_value[result$gene_site == target]-result$true_value[result$gene_site == target])^2))
        if(!target %in% names(best_q2)){
          names(new_q2) <- target
          names(new_k) <- target
          names(new_pnum) <- target
          names(new_rmse) <- target
          best_k <- c(best_k,new_k)
          best_pred_num <- c(best_pred_num,new_pnum)
          best_q2 <- c(best_q2,new_q2)
          best_rmse <- c(best_rmse,new_rmse)
        }else if(best_q2[target] < new_q2){
          best_q2[target] <- new_q2
          best_pred_num[target] <- new_pnum
          best_k[target] <- new_k
          best_rmse[target] <- new_rmse
        }
      }
    }
  }
  return(cbind(best_pred_num,best_k,best_q2,best_rmse))
}
dist_analysis2 <- function(result,cutoff,target_limit){
  m_q2 <- sapply(result,function(x)median(as.numeric(x[x[,3]>cutoff & rownames(x) %in% target_limit,3])))
  m_rmse <- sapply(result,function(x)median(as.numeric(x[x[,3]>cutoff & rownames(x) %in% target_limit,4])))
  m_npred <- sapply(result,function(x)median(as.numeric(x[x[,3]>cutoff & rownames(x) %in% target_limit,1])))
  valid_ind <- !is.na(m_q2) & !is.na(m_npred)
  return(list('q2'=m_q2[valid_ind],'rmse'=m_rmse[valid_ind],'m_pred'=m_npred[valid_ind]))
}
en_all_result <- dist_analysis(en_all_path,output_path)
dist_result <- list()
heat_result <- list()
for (i in 2:length(configs)){
  dist_result[[i]] <- dist_analysis(configs[1:i],output_path)
}
for (i in 1:length(configs_heat)){
  heat_result[[i]] <- dist_analysis(configs_heat[1:i],output_path)
}

target_limit <- Reduce(intersect,list(rownames(dist_result[[2]]),rownames(heat_result[[2]]),rownames(en_all_result)))
en_all_result_intersect <- en_all_result[target_limit,]
for(i in rownames(en_all_result_intersect)){
  en_all_result_intersect[i,3] <- max(dist_result[[length(dist_result)]][i,3],heat_result[[length(heat_result)]][i,3],en_all_result_intersect[i,3])
  en_all_result_intersect[i,4] <- min(dist_result[[length(dist_result)]][i,4],heat_result[[length(heat_result)]][i,4],en_all_result_intersect[i,4])
}


cutoff=0
dist_m <- dist_analysis2(dist_result,cutoff,target_limit)
heat_m <- dist_analysis2(heat_result,cutoff,target_limit)
#report fig 6
par(mfrow=c(1,2))
plot(log(dist_m$m_pred,base = 10),dist_m$q2,xlim = log(c(min(dist_m$m_pred,heat_m$m_pred),max(dist_m$m_pred,heat_m$m_pred)),base = 10),ylim = c(min(dist_m$q2,heat_m$q2),max(dist_m$q2,heat_m$q2)),xlab = expression('Log'[10]*' median number of predictors'),ylab = 'Median Q-squared',main = 'Q-squared vs. Predictor number')
lines(log(dist_m$m_pred,base = 10),dist_m$q2,col='blue')
points(log(heat_m$m_pred,base = 10),heat_m$q2)
lines(log(heat_m$m_pred,base = 10),heat_m$q2,col='red')
legend('topleft',legend = c('heat diffussion','distance'),col = c('red','blue'),lty = 1,bty = 'n',cex = 0.8)

plot(log(dist_m$m_pred,base = 10),dist_m$rmse,xlim = log(c(min(dist_m$m_pred,heat_m$m_pred),max(dist_m$m_pred,heat_m$m_pred)),base = 10),ylim = c(min(dist_m$rmse,heat_m$rmse),max(dist_m$rmse,heat_m$rmse)),xlab = expression('Log'[10]*' median number of predictors'),ylab = 'Median RMSE',main = 'RMSE vs. Predictor number')
lines(log(dist_m$m_pred,base = 10),dist_m$rmse,col='blue')
points(log(heat_m$m_pred,base = 10),heat_m$rmse)
lines(log(heat_m$m_pred,base = 10),heat_m$rmse,col='red')
legend(3.2,2.813,legend = c('heat diffussion','distance'),col = c('red','blue'),lty = 1,bty = 'n',cex = 0.8)

#report fig 4
mpred_dist <- c()
mpred_heat <- c()
for (i in 1:length(configs)){
  result <- result_table(paste(output_path,configs[i],sep = ''))
  mpred_dist <- c(mpred_dist,median(as.numeric(result$predictors_given[result$best_outer_q2>0])))
}
for (i in 1:length(configs_heat)){
  result <- result_table(paste(output_path,configs_heat[i],sep = ''))
  mpred_heat <- c(mpred_heat,median(as.numeric(result$predictors_given[result$best_outer_q2>0])))
}

plot(dist_range[-1],mpred_dist[-1],xlab = 'Shortest distance',ylab = 'Median number of predictors',main = 'Distance as cutoff',ylim = c(0,13000))
lines(dist_range[-1],mpred_dist[-1])
k_range <- as.numeric(paste('0.',heat_k,sep = ''))
plot(log(k_range,base = 10),mpred_heat,xlim = c(max(log(k_range,base = 10)),min(log(k_range,base = 10))),xlab = expression('Log'[10]*' log percent of heat in target node'),ylab = 'Median number of predictors',main = 'Heat diffusion as cutoff',ylim = c(0,13000))
lines(log(k_range,base = 10),mpred_heat)

