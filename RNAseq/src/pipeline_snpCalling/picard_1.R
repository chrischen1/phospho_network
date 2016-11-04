samDir  = '/groups/sorger/cchris/bam_2pass/'
outDir    = '/groups/sorger/cchris/bam_files/'
bsub_fig  = 'bsub -q short -W 4:00 -R "rusage[mem=18000]" -o '
picard_loc = '/groups/sorger/cchris/picard-tools-2.5.0/picard.jar'
  
dir.create(outDir,showWarnings = F)
sam_files <- grep('Aligned.out.sam',list.files(samDir),value = T,fixed = T)

for (i in sam_files){
  sample_id = gsub('Aligned.out.sam','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out java -jar ',picard_loc,' AddOrReplaceReadGroups I=',samDir,i,' O=',outDir,sample_id,'.bam SO=coordinate RGLB=1 RGPL=illumina RGPU=unit1 RGSM=',sample_id,sep = ''))
}