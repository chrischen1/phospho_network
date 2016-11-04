xxx=read.delim("~/Documents/workspace/phospho_network/RNAseq/annotation/20160325mRNAseq1key.tsv", header=FALSE, comment.char="#", stringsAsFactors=FALSE,skip = 2)

yyy=read.delim("~/Documents/workspace/phospho_network/RNAseq/annotation/20160404mRNAseq2key.tsv", header=FALSE, comment.char="#", stringsAsFactors=FALSE)

all <- rbind(xxx,yyy[-1,])

all2 <- apply(all,2,function(x)gsub('\xa0','',x))
rep_tag <- gsub('[^a-z+]','',all2[,2])
sample_id <- paste(gsub('\\D','',all2[,2]),'_S',all2[,1],sep = '')
cellines <- toupper(gsub(' ','',all2[,4]))
cellines <- gsub('MDA','MDAMB',cellines)
cellines_annotation <- cbind(sample_id,cellines,rep_tag)
write.table(cellines_annotation,'~/Documents/workspace/phospho_network/RNAseq/annotation/RNAseq_sample.tsv',sep = '\t',quote = F,row.names = F)
