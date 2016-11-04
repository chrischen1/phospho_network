bamDir  = '/groups/sorger/cchris/bam_files_dupSplit_realign/'
outDir    = '/groups/sorger/cchris/bam_files_dupSplit_realign_BQSR/'
bsub_fig  = 'bsub -q short -W 12:00 -R "rusage[mem=28000]" -o '
gatk_loc = '/groups/sorger/cchris/gatk/GenomeAnalysisTK.jar'
ref_fasta = '/groups/sorger/cchris/external_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

dir.create(outDir,showWarnings = F)
bam_files <- grep('bam',list.files(bamDir),value = T,fixed = T)

for (i in bam_files){
  sample_id = gsub('.bam','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out java -jar ',gatk_loc,' -T PrintReads -R ',ref_fasta,' -I ',bamDir,i,' -o ',outDir,i,' -BQSR ',outDir,sample_id,'.grp',sep = ''))
}
