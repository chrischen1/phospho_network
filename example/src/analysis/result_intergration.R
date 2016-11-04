result_path = '~/Documents/workspace/phospho_network/example/script_files/output/'
output_path = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/'
k_range = c('0','1','01','001','0008','0006','0004','0002','0001','00001','000001')
config_prefix = 'parallel_config_en_large_heat_ud_k'
config_prefix_rf = 'rf/rfL_heat_ud_results_k'
cal_q2 <- function(x,y){
  return(round(1-sum((x-y)^2)/sum((y-mean(y))^2),4))
}
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

result_integration <- function(result_path,k_range,config_prefix){
  result_best <- NULL
  beta_best <- NULL
  for(k in k_range){
    k_val <- as.numeric(paste('0.',k,sep = ''))
    config_path <- paste(result_path,paste(config_prefix,k,sep = ''),sep = '')
    result <- get_result(config_path,keyword = 'result')
    beta <- get_result(config_path,keyword = 'beta')
    targets <- unique(result$gene_site)
    for(target in targets){
      if(is.null(result_best)){
        result_best <- cbind(result[result$gene_site==target,],'k'=k_val)
        beta_best <- beta[beta$gene_site==target,]
        }
      else if(!target %in% unique(result_best$gene_site)){
        result_best <- rbind(result_best,cbind(result[result$gene_site == target,],'k'=k_val))
        beta_best <- rbind(beta_best,beta[beta$gene_site == target,])
      }
      else if (unique(result$best_outer_q2[result$gene_site==target]) > unique(result_best$best_outer_q2[result_best$gene_site==target])){
        result_best <- rbind(result_best[result_best$gene_site != target,],cbind(result[result$gene_site == target,],'k'=k_val))
        beta_best <- rbind(beta_best[beta_best$gene_site != target,],beta[beta$gene_site == target,])
      }
    }
  }
  return(list('result_best'=result_best,'beta_best'=beta_best))
}

integration_heat <- result_integration(result_path,k_range[2:10],config_prefix)
integration_rf <- result_integration(result_path,k_range[2:10],config_prefix_rf)
write.csv(integration_heat$beta_best,paste(output_path,'beta_heat_best_en_00001.csv',sep = ''),row.names = F)
write.csv(integration_heat$result_best,paste(output_path,'result_heat_best_en_00001.csv',sep = ''),row.names = F)
write.csv(integration_rf$beta_best,paste(output_path,'beta_heat_best_rf_00001.csv',sep = ''),row.names = F)
write.csv(integration_rf$result_best,paste(output_path,'result_heat_best_rf_00001.csv',sep = ''),row.names = F)

integration_heat2 <- result_integration(result_path,k_range[1:10],config_prefix)
integration_rf2 <- result_integration(result_path,k_range[1:10],config_prefix_rf)
write.csv(integration_heat2$beta_best,paste(output_path,'beta_heat_best_en_toall.csv',sep = ''),row.names = F)
write.csv(integration_heat2$result_best,paste(output_path,'result_heat_best_en_toall.csv',sep = ''),row.names = F)
write.csv(integration_rf2$beta_best,paste(output_path,'beta_heat_best_rf_toall.csv',sep = ''),row.names = F)
write.csv(integration_rf2$result_best,paste(output_path,'result_heat_best_rf_toall.csv',sep = ''),row.names = F)
