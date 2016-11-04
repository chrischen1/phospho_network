rpkm_file1  = "~/Documents/workspace/phospho_network/RNAseq/processed_data/star/rpkm.csv"
rpkm_file2  = "~/Documents/workspace/phospho_network/RNAseq/processed_data/test.tsv"
rpkm_file1_out  = "~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_exp.tsv"
rpkm_file2_out  = "~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_paper3080.tsv"

library(org.Hs.eg.db)
library(annotate)
rpkm1 <- read.csv(rpkm_file1, stringsAsFactors=FALSE)
rpkm1 <- rpkm1[!rpkm1$X %in% rpkm1$X[duplicated(rpkm1$X)],]
rownames(rpkm1) <- rpkm1$X
rpkm1 <- rpkm1[,-c(1,40,41)]
colnames(rpkm1) <- c('BT-549','HCC38','HCC70','MDA-MB-134-VI','MDA-MB-157','MDA-MB-361','MDA-MB-436','MDA-MB-453','MDA-MB-468','PDX_HCI002',
                     'PDX_HCI002_techRep','HCC1806','HCC1143','HCC1395','HCC1419','HCC1500','HCC1937','HCC1954','CAL-120','CAL-51','CAL-85','SUM149PT',
                     'MCF 10A_bioRep','MCF 10A','SUM159PT','SUM1315MO2','T47D','CAMA-1','HME1','HCC1428','BT-20','SK-BR-3','Hs 578T','MDA-MB-231',
                     'MCF7','PDX1258','PDX1258_techRep','PDX1328','PDX1328_techRep')
write.table(rpkm1,rpkm_file1_out,sep = '\t',row.names = T,col.names = T)

rpkm2 <- read.delim(rpkm_file2,stringsAsFactors=FALSE)
cell_info <- read.delim("~/Documents/workspace/phospho_network/RNAseq/cellines_info/celllines_info_Genentech.txt",stringsAsFactors = F)
breast_cells <- cell_info$Cell.line[cell_info$Primary.Tissue == 'Breast']
rpkm2_genes <- getSYMBOL(as.character(rpkm2$GeneID), data='org.Hs.eg')
rpkm2$GeneID <- rpkm2_genes
rpkm2 <- rpkm2[!rpkm2$GeneID %in% rpkm2$GeneID[duplicated(rpkm2$GeneID)],]
rownames(rpkm2) <- rpkm2$GeneID
rpkm2 <- rpkm2[,colnames(rpkm2) %in% breast_cells]
colnames(rpkm2) <- c('CAL-51','DU4475','HCC2157','HCC2688','HCC2911','HCC1493','MCF10DCIS.com','HCC1143','HCC1187',
                     'HCC1395','HCC1419','HCC1428','HCC1500','HCC1569','HCC1599','HCC1806','HCC1937','HCC1954','HCC202',
                     'HCC2218','HCC38','HCC70')
write.table(rpkm2,rpkm_file2_out,sep = '\t',row.names = T,col.names = T)








