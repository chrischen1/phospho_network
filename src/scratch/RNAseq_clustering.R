pdf('with_subtype.pdf')
rpkm <- read.delim("~/Downloads/rpkm.tsv", row.names=1)
cell_line_categories <- read.delim("~/Downloads/cell_line_categories.tsv",as.is = T)

cell_type = paste(cell_line_categories$Receptor.Status,cell_line_categories$Molecular.Subtype)
cell_type[cell_line_categories$Molecular.Subtype==''] <- paste(cell_type[cell_line_categories$Molecular.Subtype==''],cell_line_categories$Cell.Line[cell_line_categories$Molecular.Subtype==''])
cell_type2<-gsub(' +','_',cell_type)
names(cell_type2)=cell_line_categories$Cell.Line

na_cell =colnames(rpkm)[!colnames(rpkm) %in% names(cell_type2)]
na_type <- rep('NM_Non_malignant,_Basal',length(na_cell))
names(na_type) <- na_cell
cell_type2 <- c(cell_type2, na_type)

x=rpkm[apply(rpkm,1,var)>0,]
logx=log(x+0.0001)
colnames(logx)<-paste(cell_type2[colnames(logx)],colnames(logx),sep = '_')
rpkm_all_cor <- cor(logx,method = 'spearman')
dissimilarity <- 1 - rpkm_all_cor
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="clustering with genes")


library(ggfortify)
pcx=prcomp(t(logx), center = TRUE, scale. = TRUE)
plot(pcx,type='l')
exp_type=cell_type2[colnames(rpkm)]


y=data.frame(cbind(pcx$x,'type'=factor(exp_type)))
aaa=cbind('cell'=colnames(rpkm),'type'=exp_type[colnames(rpkm)])
autoplot(pcx,data=aaa,colour='type',size=6)
rpkm_all_cor <- cor(t(pcx$x[,1:20]),method = 'spearman')
dissimilarity <- 1 - rpkm_all_cor
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="clustering with first 5 PCs")
dev.off()

pdf('without_subtype.pdf')
cell_type = paste(cell_line_categories$Receptor.Status)
cell_type2<-gsub(' +','_',cell_type)
names(cell_type2)=cell_line_categories$Cell.Line
cell_type2[cell_line_categories$Molecular.Subtype==''] <-cell_line_categories$Cell.Line[cell_line_categories$Molecular.Subtype=='']
na_cell =colnames(rpkm)[!colnames(rpkm) %in% names(cell_type2)]
na_type <- rep('NM',length(na_cell))
names(na_type) <- na_cell
cell_type2 <- c(cell_type2, na_type)

x=rpkm[apply(rpkm,1,var)>0,]
logx=log(x+0.0001)
colnames(logx)<-paste(cell_type2[colnames(logx)],colnames(logx),sep = '_')
rpkm_all_cor <- cor(logx,method = 'spearman')
dissimilarity <- 1 - rpkm_all_cor
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="clustering with genes")


library(ggfortify)
pcx=prcomp(t(logx), center = TRUE, scale. = TRUE)
plot(pcx,type='l')
exp_type=cell_type2[colnames(rpkm)]


y=data.frame(cbind(pcx$x,'type'=factor(exp_type)))
aaa=cbind('cell'=colnames(rpkm),'type'=exp_type[colnames(rpkm)])
autoplot(pcx,data=aaa,colour='type',size=6)
rpkm_all_cor <- cor(t(pcx$x[,1:20]),method = 'spearman')
dissimilarity <- 1 - rpkm_all_cor
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="clustering with first 5 PCs")
dev.off()
