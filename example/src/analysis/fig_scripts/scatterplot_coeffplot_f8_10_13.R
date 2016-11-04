result_file_heat = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_en_00001.csv'
result_file_all = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_all.csv'
result_file_rf = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_rf_00001.csv'
beta_file_heat   = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_heat_best_en_00001.csv'
beta_file_all   = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_all.csv'
beta_file_en   = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_rf_00001.csv'
beta_file_rf   = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_heat_best_rf_00001.csv'
nbeta = 10
id    = 'EGFR_1068'

result_heat <- read.csv(result_file_heat,as.is = T)
beta_heat   <- read.csv(beta_file_heat,as.is = T)
result_all <- read.csv(result_file_all,as.is = T)
beta_all   <- read.csv(beta_file_all,as.is = T)

result_rf <- read.csv(result_file_rf,as.is = T)
beta_rf   <- read.csv(beta_file_rf,as.is = T)


plot_result <- function(unique_id,result,beta,nbeta = 10,note = ''){
  predict_matrix <- result[result$gene_site == unique_id,]
  beta_matrix    <- beta[beta$gene_site == unique_id,]
  true_val <- result[result$gene_site == unique_id,'true_value']
  pred_val <- result[result$gene_site == unique_id,'predict_value']
  q2_value <- unique(result[result$gene_site == unique_id,'best_outer_q2'])
  # pval <- round(cor.test(pred_val,true_val,alternative = 'greater')$p.value,5)
  # if(pval < 1e-6){pval <- '<=1e-6'}else{pval <- round(pval,5)}
  smoothScatter(true_val,pred_val,ylab = paste('Predicted,q-sqad =',round(q2_value,4)),xlab = 'True',main = paste('Regression for ',unique_id,note))
  beta_matrix2 <- t(beta_matrix[,-(1:2)])
  colnames(beta_matrix2) <- beta_matrix$predictor
  plot_names <- names(tail(sort(apply(beta_matrix2,2,function(x)abs(median(x)))),nbeta))
  boxplot(beta_matrix2[,names(sort(apply(beta_matrix2[,plot_names],2,median)))],cex.axis=0.45,las = 2, main = paste('coefficients in CV',note),
          ylab = paste('predictors used:',ncol(beta_matrix2),'/',unique(predict_matrix$predictors_given,sep = '')))
}

plot_compare <- function(result1,result2,beta1,beta2,note1='',note2='',num_plot = 1){
  q2_1 <- c()
  q2_2 <- c()
  id_intersect  <- intersect(result1$gene_site,result2$gene_site)
  for(i in id_intersect){
    q2_1 <- c(q2_1,unique(result1$best_outer_q2[result1$gene_site == i]))
    q2_2 <- c(q2_2,unique(result2$best_outer_q2[result2$gene_site == i]))
  }
  q2_diff <- q2_1-q2_2
  names(q2_diff) <- id_intersect
  q2_diff <- sort(q2_diff)
  diff_ids <- head(names(q2_diff),num_plot)
  for(id in diff_ids){
    plot_result(id,result1,beta1,note = note1)
    plot_result(id,result2,beta2,note = note2)
  }
}

#fig8a 9a
par(mfrow = c(1,2))
plot_result(id,result_heat,beta_heat)
#fig10 heat vs all 
aaa=result_heat[,c(1,6)]
aaa=unique(aaa)
aq=aaa[,2]
names(aq)=aaa[,1]

bbb=result_all[,c(1,6)]
bbb=unique(bbb)
bq=bbb[,2]
names(bq)=bbb[,1]

diffq=aq-bq[names(aq)]

#fig8a
par(mfrow = c(2,2))
plot_compare(result_all,result_heat[result_heat[,1]==names(sort(diffq,decreasing = T))[6],],beta_all,beta_heat,'All measurements','Heat diffusion')
plot_compare(result_heat,result_all,beta_heat,beta_all,'Heat diffusion','All measurements')
#fig13 heat vs rf
plot_compare(result_heat[result_heat[,1]=='SRC_416',],result_rf,beta_heat,beta_rf,'Elastic Net','Random Forest')


