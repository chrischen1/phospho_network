#first merge the replicates by average, and perform log2(x+1) transform
rpkm_raw_file = '~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_lincs.tsv'
out_file      = '~/Documents/workspace/phospho_network/RNAseq/processed_data/rpkm_processed.tsv'

rpkm_raw <- read.delim(rpkm_raw_file,as.is = T)[-1,]
genes <- rownames(rpkm_raw)
rpkm_raw <- apply(rpkm_raw,2,as.numeric)
rownames(rpkm_raw) <- genes
colnames(rpkm_raw) <- gsub('\\.\\d+','',colnames(rpkm_raw))
samples <- unique(colnames(rpkm_raw))

rpkm <- NULL
for (sid in samples){
  if(sum(sid==colnames(rpkm_raw))==1){
    rpkm <- cbind(rpkm,rpkm_raw[,sid])
  }else{
    avg_result <- apply(rpkm_raw[,colnames(rpkm_raw) == sid],1,mean)
    rpkm <- cbind(rpkm,avg_result)
  }
}
rpkm <- log(rpkm+1,base = 2)
colnames(rpkm) <- samples
rownames(rpkm) <- genes
write.table.hdr <- function( X, file, idtxt="id" )
{
  cat( idtxt, "\t", sep="", file=file )
  suppressWarnings(utils::write.table( X, file=file, quote=FALSE, sep="\t", append=TRUE )) 
}


write.table.hdr(rpkm,out_file)
