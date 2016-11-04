bamDir  = '/groups/sorger/cchris/bam_files_markdup/'
outDir    = '/groups/sorger/cchris/bam_files_dupSplit/'
bsub_fig  = 'bsub -q short -W 12:00 -R "rusage[mem=18000]" -o '
gatk_loc = '/groups/sorger/cchris/gatk/GenomeAnalysisTK.jar'
ref_fasta = '/groups/sorger/cchris/external_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

dir.create(outDir,showWarnings = F)
bam_files <- grep('bam',list.files(bamDir),value = T,fixed = T)

for (i in bam_files){
  sample_id = gsub('.bam','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out java -jar ',gatk_loc,' -T SplitNCigarReads -R ',ref_fasta,' -I ',bamDir,i,' -o ',outDir,i,' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS',sep = ''))
}