rpkm_file1  = "~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_exp.tsv"
rpkm_file2  = "~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_paper3080.tsv"

library(mvtnorm)
rpkm1 <- read.delim(rpkm_file1,stringsAsFactors = F)
rpkm2 <- read.delim(rpkm_file2,stringsAsFactors = F)

var_genes_rpkm1 <- rownames(rpkm1)[apply(rpkm1,1,var)>0]
var_genes_rpkm2 <- rownames(rpkm2)[apply(rpkm2,1,var)>0]
var_genes_intersect <- intersect(var_genes_rpkm1,var_genes_rpkm2)

rpkm1i <- rpkm1[var_genes_intersect,]
rpkm2i <- rpkm2[var_genes_intersect,]
colnames(rpkm2i) <- paste(colnames(rpkm2i),'ref',sep = '_')

#spearman correlation
rpkm_all <- cbind(rpkm1i,rpkm2i)
rpkm_all_cor <- cor(rpkm_all,method = 'spearman')
dissimilarity <- 1 - rpkm_all_cor
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="")

#sva
source('~/Documents/workspace/phospho_network/RNAseq/src/combat2.R')
rpkm_norm <- combat2(log(rpkm1i+1,base = 2),log(rpkm2i+1,base = 2))
rpkm_all_norm <- cbind(rpkm_norm[[1]],rpkm_norm[[2]])
var_1kgenes <- names(tail(sort(apply(rpkm_all_norm,1,var)),1000))
heatmap_m <- rpkm_all_norm[var_1kgenes,]
row_ind <- order.dendrogram(as.dendrogram(hclust(dist(heatmap_m))))
heatmap(as.matrix(heatmap_m[row_ind,]),Rowv = NA,labRow = NA)

write.table(rpkm1[var_genes_rpkm1,],"~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_exp_non0.tsv",sep = '\t',row.names = T,col.names = T)
write.table(rpkm2[var_genes_rpkm2,],"~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_paper3080_non0.tsv",sep = '\t',row.names = T,col.names = T)

tpm_exp_non0 <- apply(rpkm1[var_genes_rpkm1,],2,function(x)10^6*x/sum(x))
tpm_paper3080_non0 <- apply(rpkm2[var_genes_rpkm2,],2,function(x)10^6*x/sum(x))
write.table(tpm_exp_non0,"~/Documents/workspace/phospho_network/RNAseq/processed_data/tpm_exp_non0.tsv",sep = '\t',row.names = T,col.names = T)
write.table(tpm_paper3080_non0,"~/Documents/workspace/phospho_network/RNAseq/processed_data/tpm_paper3080_non0.tsv",sep = '\t',row.names = T,col.names = T)

rpkms_1kgenes <- rpkm_all[var_1kgenes,]
rpkms_1kgenes_norm <- rpkm_all_norm[var_1kgenes,]
tpm_1kgenes <- cbind(tpm_exp_non0[var_1kgenes,],tpm_paper3080_non0[var_1kgenes,])

smoothScatter(log(as.numeric(as.matrix(rpkms_1kgenes_norm))),log(as.numeric(as.matrix(rpkms_1kgenes))),xlab = 'RPKM Normalized',ylab = 'RPKM')
smoothScatter(log(as.numeric(as.matrix(rpkms_1kgenes_norm))),log(as.numeric(as.matrix(tpm_1kgenes))),xlab = 'RPKM Normalized',ylab = 'TPM')
smoothScatter(log(as.numeric(as.matrix(tpm_1kgenes))),log(as.numeric(as.matrix(rpkms_1kgenes))),xlab = 'TPM',ylab = 'RPKM')

par(mfrow = c(3,1))
plot(density(as.numeric(as.matrix(rpkm_all))),main = 'Density for orginal RPKM')
plot(density(as.numeric(as.matrix(log(rpkm_all+1,2)))),main = 'Density for log2 RPKM')
plot(density(as.numeric(as.matrix(rpkm_all_norm))),main = 'Density for Normalized log2 RPKM')


