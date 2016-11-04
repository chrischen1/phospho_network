library(ggfortify)
rpkm <- read.delim("~/Downloads/rpkm.tsv", row.names=1)
cell_line_categories <- read.delim("~/Downloads/cell_line_categories.tsv",as.is = T)

pc_plot <- function(logx,cell_line,t,cells){
  pcx=prcomp(t(logx), center = TRUE, scale. = TRUE)
  g<-autoplot(pcx,data=cell_line,colour=t,size=6)
  pc_sum <- summary(pcx)$importance
  g=g+geom_text(aes(label=cells), size=1)+labs(x = paste("PC1",round(pc_sum[2,1],4)),y=paste('PC2',round(pc_sum[2,2],4)))
  g
}

cell_line_categories <- rbind(cell_line_categories,c('HME1','NM','Non malignant, Basal'))
rownames(cell_line_categories) <- cell_line_categories$Cell.Line
cell_line <- cell_line_categories[colnames(rpkm),]
cell_line <- cbind(cell_line,'type'=gsub(' ','_',paste(cell_line$Receptor.Status,cell_line$Molecular.Subtype)))
logx <- log(rpkm[apply(rpkm,1,var)>0,]+0.0001,base = 2)
cells <- gsub('.+_.+_(.+)','\\1',colnames(logx))

# rpkm_all_cor <- cor(logx,method = 'spearman')
# dissimilarity <- 1 - rpkm_all_cor
# distance <- as.dist(dissimilarity)
# plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="clustering with genes")

pdf('with_all_subtype.pdf')
pc_plot(logx,cell_line,'type',cells)
dev.off()

pdf('with_Receptor_Status.pdf')
pc_plot(logx,cell_line,'Receptor.Status',cells)
dev.off()

pdf('with_Molecular_Subtype.pdf')
pc_plot(logx,cell_line,'Molecular.Subtype',cells)
dev.off()

