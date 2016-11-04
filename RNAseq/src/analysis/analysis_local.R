library(corrplot)
library(edgeR)


subread_path = "~/Documents/workspace/phospho_network/RNAseq/processed_data/subread/"
star_path = "~/Documents/workspace/phospho_network/RNAseq/processed_data/star/"

rep_matrix = matrix(c("19","19_tech_replic","8","8_tech_replic","9","9_tech_replic","2","2_bio_replic"),ncol = 2,byrow = T)

build_deg <- function(result_path){
  counts <- read.csv(paste(result_path,"count_table.csv",sep = ''),row.names=1)
  annotation <-read.csv(paste(result_path,"annotation_table.csv",sep = ''),row.names=1,as.is = T)
  stat <- read.csv(paste(result_path,"stats_table.csv",sep = ''),row.names=1)
  targets <- gsub('.','_',colnames(counts),fixed = T)
  targets <- gsub('.+files_(\\w+)_S.+','\\1',targets)
  colnames(counts) <- targets
  colnames(stat)[-1] <- targets
  fc <- list('counts'=counts,'annotation'=annotation,'targets'=targets,'stat'=stat)
  return(fc)
}

par(mfrow = c(2,2))

fc <- build_deg(subread_path)
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
x_rpkm <- rpkm(x,x$genes$Length)
cor_m <- cor(fc$counts[,as.character(rep_matrix)])
corrplot(cor_m,type="upper",title = 'Subread raw count')
cor_m2 <- cor(x_rpkm[,as.character(rep_matrix)])
corrplot(cor_m2,type="upper",title = 'Subread RPKM')

fc2 <- build_deg(star_path)
x <- DGEList(counts=fc2$counts, genes=fc2$annotation[,c("GeneID","Length")])
x_rpkm <- rpkm(x,x$genes$Length)
cor_m <- cor(fc2$counts[,as.character(rep_matrix)])
corrplot(cor_m,type="upper",title = 'STAR raw count')
cor_m2 <- cor(x_rpkm[,as.character(rep_matrix)])
corrplot(cor_m2,type="upper",title = 'STAR RPKM')


for(i in 1:nrow(rep_matrix)){
  plot(x_rpkm[,rep_matrix[i,1]],x_rpkm[,rep_matrix[i,2]],xlab = rep_matrix[i,1],ylab = rep_matrix[i,2],main = 'RPKM raw')
}

for(i in 1:nrow(rep_matrix)){
  plot(log(x_rpkm[,rep_matrix[i,1]]+1),log(x_rpkm[,rep_matrix[i,2]]+1),xlab = rep_matrix[i,1],ylab = rep_matrix[i,2],main = 'RPKM log')
}
for(i in 1:nrow(rep_matrix)){
  plot(fc2$counts[,rep_matrix[i,1]],fc2$counts[,rep_matrix[i,2]],xlab = rep_matrix[i,1],ylab = rep_matrix[i,2],main = 'counts raw')
}
for(i in 1:nrow(rep_matrix)){
  plot(log(fc2$counts[,rep_matrix[i,1]]+1),log(fc2$counts[,rep_matrix[i,2]]+1),xlab = rep_matrix[i,1],ylab = rep_matrix[i,2],main = 'counts log')
}

