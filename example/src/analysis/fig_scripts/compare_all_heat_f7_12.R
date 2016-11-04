en_heat_best_file = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_en_00001.csv'
rf_heat_best_file = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_rf_00001.csv'
en_heat_best_beta = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_heat_best_en_00001.csv'
rf_heat_best_beta = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_heat_best_rf_00001.csv'
all_file          = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_all.csv'
all_beta_file     = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/beta_all.csv'
require(gridExtra)
require(ggplot2)
cal_q2 <- function(x,y){
  return(round(1-sum((x-y)^2)/sum((y-mean(y))^2),4))
}

result_heat_en <- read.csv(en_heat_best_file,as.is = T)
result_all <- read.csv(all_file,as.is = T)
beta_en <- read.csv(en_heat_best_beta,as.is = T)
beta_all <- read.csv(all_beta_file,as.is = T)
# fig7
target_intersect <- Reduce(intersect,list(result_heat_en$gene_site,result_all$gene_site))
result_heat_en <- result_heat_en[result_heat_en$gene_site %in%target_intersect,]
result_all <- result_all[result_all$gene_site %in%target_intersect,]
par(mfrow = c(1,2))
smoothScatter(result_all$true_value,result_all$predict_value,xlab = 'True',xlim = c(0,1),ylim = c(0,1),
              ylab = paste('Predicted, q-sqad =',cal_q2(result_all$predict_value,result_all$true_value)),main = 'Regression based on all predictors')
smoothScatter(result_heat_en$true_value,result_heat_en$predict_value,xlab = 'True',xlim = c(0,1),ylim = c(0,1),
              ylab = paste('Predicted, q-sqad =',cal_q2(result_heat_en$predict_value,result_heat_en$true_value)),main = 'Regression based on heat diffusion')

df_q2 <- data.frame(rbind(cbind(result_all$predictors_given,result_all$best_outer_q2,'All measurements'),
                          cbind(result_heat_en$predictors_given,result_heat_en$best_outer_q2,'Heat diffusion')))
colnames(df_q2) <- c('Predictor_Number','Q_squared','pred_selection')
df_q2$Q_squared <- as.numeric(as.character(df_q2$Q_squared))
df_q2$Predictor_Number <- as.numeric(as.character(df_q2$Predictor_Number))
df_q2_heat <- df_q2[df_q2$pred_selection == 'Heat diffusion',]
p <- ggplot(df_q2,aes(pred_selection,Q_squared))+ geom_boxplot()+labs(title = "Prediction with different predictor selection")
pred_use_en <- cbind(beta_en[,1],apply(beta_en[,-(1:2)],1,median))
pred_use_en_table <- table(pred_use_en[pred_use_en[,2]!=0,1])
pred_use_all <- cbind(beta_all[,1],apply(beta_all[,-(1:2)],1,median))
pred_use_all_table <- table(pred_use_all[pred_use_all[,2]!=0,1])
df_npred <- data.frame(rbind(cbind(as.numeric(pred_use_en_table),'scope'='Heat diffusion'),cbind(as.numeric(pred_use_all_table),'scope'='All measurements')))
df_npred$V1 <- log(as.numeric(as.character(df_npred$V1)),2)
names(df_npred)[1] <- 'Predictor_Number'
p2 <- ggplot(df_npred, aes(Predictor_Number,fill = scope,color=scope)) +geom_density(alpha = 0.1)+labs(title = "Number of predictors selected for prediction")
grid.arrange(p, p2, ncol=2)
all_q <-unique(df_q2$Q_squared[df_q2$pred_selection=='All measurements'])
heat_q <- unique(df_q2$Q_squared[df_q2$pred_selection=='Heat diffusion'])
wilcox.test(all_q,heat_q)
all_npred <- unique(df_npred$Predictor_Number[df_npred$scope=='All measurements'])
heat_npred <- unique(df_npred$Predictor_Number[df_npred$scope=='Heat diffusion'])
wilcox.test(all_npred,heat_npred)
#fig12 top
result_heat_en <- read.csv(en_heat_best_file,as.is = T)
result_heat_rf <- read.csv(rf_heat_best_file,as.is = T)
target_intersect <- Reduce(intersect,list(result_heat_en$gene_site,result_heat_rf$gene_site))
result_heat_en <- result_heat_en[result_heat_en$gene_site %in%target_intersect,]
result_heat_rf <- result_heat_rf[result_heat_rf$gene_site %in%target_intersect,]
smoothScatter(result_heat_rf$true_value,result_heat_rf$predict_value,xlab = 'True',xlim = c(0,1),ylim = c(0,1),
              ylab = paste('Predicted, q-sqad =',cal_q2(result_heat_rf$predict_value,result_heat_rf$true_value)),main = 'Regression based on random forest')

df_q2 <- data.frame(rbind(cbind(result_heat_rf$predictors_given,result_heat_rf$best_outer_q2,'random forest'),
                          cbind(result_heat_en$predictors_given,result_heat_en$best_outer_q2,'elastic net')))
colnames(df_q2) <- c('Predictor_Number','Q_squared','pred_selection')
df_q2$Q_squared <- as.numeric(as.character(df_q2$Q_squared))
df_q2$Predictor_Number <- as.numeric(as.character(df_q2$Predictor_Number))
df_q2_en <- df_q2[df_q2$pred_selection == 'elastic net',]
df_q2_rf <- df_q2[df_q2$pred_selection == 'random forest',]
p <- ggplot(df_q2,aes(pred_selection,Q_squared))+ geom_boxplot()+labs(title = "Prediction with different algorithms")
#fig12 bottom
beta_rf <- read.csv(rf_heat_best_beta,as.is = T)
pred_use_rf <- cbind(beta_rf[,1],apply(beta_rf[,-(1:2)],1,median))
pred_use_rf_table <- table(pred_use_rf[pred_use_rf[,2]!=0,1])
df_npred <- data.frame(rbind(cbind(as.numeric(pred_use_en_table),'algos'='elastic net'),cbind(as.numeric(pred_use_rf_table),'algos'='random forest')))
df_npred$V1 <- as.numeric(as.character(df_npred$V1))
names(df_npred)[1] <- 'Predictor_Number'
p2 <- ggplot(df_q2, aes(Predictor_Number,fill = pred_selection,color=pred_selection)) +geom_density(alpha = 0.1)+labs(title = "Number of input for training")
p3 <- ggplot(df_npred, aes(Predictor_Number,fill = algos,color=algos)) +geom_density(alpha = 0.1)+labs(title = "Number of predictors selected for prediction")
grid.arrange(p2, p3, ncol=2)
rf_q <-unique(df_q2$Q_squared[df_q2$pred_selection=='random forest'])
en_q <- unique(df_q2$Q_squared[df_q2$pred_selection=='elastic net'])
rf_predt <-unique(df_q2$Predictor_Number[df_q2$pred_selection=='random forest'])
en_predt <- unique(df_q2$Predictor_Number[df_q2$pred_selection=='elastic net'])
rf_npred <- unique(df_npred$Predictor_Number[df_npred$algos=='random forest'])
en_npred <- unique(df_npred$Predictor_Number[df_npred$algos=='elastic net'])
wilcox.test(rf_q,en_q)
wilcox.test(rf_npred,en_npred)
wilcox.test(rf_predt,en_predt)



