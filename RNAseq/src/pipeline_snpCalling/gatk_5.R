bamDir  = '/groups/sorger/cchris/bam_files_dupSplit_realign_BQSR/'
outDir    = '/groups/sorger/cchris/vcf_files/'
bsub_fig  = 'bsub -q long -W 48:00 -R "rusage[mem=18000]" -o '
gatk_loc = '/groups/sorger/cchris/gatk/GenomeAnalysisTK.jar'
ref_fasta = '/groups/sorger/cchris/external_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

dir.create(outDir,showWarnings = F)
bam_files <- grep('bam',list.files(bamDir),value = T,fixed = T)

for (i in bam_files){
  sample_id = gsub('.bam','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out java -jar ',gatk_loc,' -T HaplotypeCaller -R ',ref_fasta,' -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -I ',bamDir,i,' -o ',outDir,sample_id,'.vcf',sep = ''))
}