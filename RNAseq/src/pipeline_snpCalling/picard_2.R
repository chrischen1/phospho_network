bamDir  = '/groups/sorger/cchris/bam_files/'
outDir    = '/groups/sorger/cchris/bam_files_markdup/'
bsub_fig  = 'bsub -q short -W 4:00 -R "rusage[mem=48000]" -o '
picard_loc = '/groups/sorger/cchris/picard-tools-2.5.0/picard.jar'

dir.create(outDir,showWarnings = F)
bam_files <- grep('bam',list.files(bamDir),value = T,fixed = T)

for (i in bam_files){
  sample_id = gsub('.bam','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out java -jar ',picard_loc,' MarkDuplicates I=',bamDir,i,' O=',outDir,sample_id,'.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=',outDir,sample_id,'.metrics',sep = ''))
}