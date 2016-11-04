temp_path = '~/Documents/workspace/phospho_network/example/script_files/output'
rna_result_file = '~/Documents/workspace/phospho_network/example/script_files/output/pred_rna.csv'
config_path = "~/Documents/workspace/phospho_network/example/report/scatterplot"
render_config = 'enL_all_results'
nbeta = 10  #max number of predictors showed in boxplot
result_path = paste(temp_path,render_config,sep = '/')

cal_q2 <- function(x,y){
  return(round(1-sum((x-y)^2)/sum((y-mean(y))^2),4))
}

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
rna_result <- read.csv(rna_result_file,as.is = T)

target_intersect <- intersect(result$gene_site,rna_result$target)
result_intersect <- result[result$gene_site %in% target_intersect,]
rna_result_intersect <- rna_result[rna_result$target %in% target_intersect,]

par(mfrow = c(1,2))
q2_rna_list <- c()
q2_list <- unique(result_intersect$best_outer_q2)



rna_target<- unique(rna_result_intersect$target)
q2_all <- c()
for (rna in rna_target){
  result_rna <- rna_result_intersect[rna_result_intersect$target == rna,]
  q2_rna_list <- c(q2_rna_list,cal_q2(result_rna$pred,result_rna$true))
  q2_all <- c(q2_all,rep(cal_q2(result_rna$pred,result_rna$true),nrow(result_rna)))
  
}
rna_result_intersect <- cbind(rna_result_intersect,'q2' = q2_all)

q2 <- cal_q2(result_intersect$predict_value,result_intersect$true_value)
q2_rna <- cal_q2(rna_result_intersect$pred,rna_result_intersect$true)

wilcox.test(unique(result_intersect$best_outer_q2[result_intersect$best_outer_q2>0]),unique(rna_result_intersect$q2[rna_result_intersect$q2>0]))

smoothScatter(result$true_value,result$predict_value,xlab = 'True',ylab = paste('Predicted, q-sqad =',cal_q2(result$predict_value,result$true_value)),main = 'Regression based on all features')
smoothScatter(rna_result$true,rna_result$pred,xlab = 'True',ylab = paste('Predicted, q-sqad =',cal_q2(rna_result$pred,rna_result$true)),main = 'Regression based on RNA only')


