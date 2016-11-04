result_file_heat = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_en_00001.csv'
result_file_all = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_all.csv'
result_file_rf = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_rf_00001.csv'
id    = 'EGFR_1068'

result_heat <- read.csv(result_file_heat,as.is = T)
result_all <- read.csv(result_file_all,as.is = T)



plot_result <- function(unique_id,result,note = ''){
  predict_matrix <- result[result$gene_site == unique_id,]
  true_val <- result[result$gene_site == unique_id,'true_value']
  pred_val <- result[result$gene_site == unique_id,'predict_value']
  q2_value <- unique(result[result$gene_site == unique_id,'best_outer_q2'])
  # pval <- round(cor.test(pred_val,true_val,alternative = 'greater')$p.value,5)
  # if(pval < 1e-6){pval <- '<=1e-6'}else{pval <- round(pval,5)}
  smoothScatter(true_val,pred_val,ylab = paste('Predicted,q-sqad =',round(q2_value,4)),xlab = 'True',main = paste('Regression for ',unique_id,note))
}

plot_compare <- function(result1,result2,note1='',note2='',num_plot = 1){
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
    plot_result(id,result1,note = note1)
    plot_result(id,result2,note = note2)
  }
}


#fig8a
par(mfrow = c(1,2))
pdf('~/egfr_scatter.pdf')
plot_compare(result_all,result_heat[result_heat$gene_site==id,],'All measurements','Heat diffusion')
dev.off()