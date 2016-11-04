vcfDir  = '/groups/sorger/cchris/vcf_files/'
outDir    = '/groups/sorger/cchris/vcf_files_filtered/'
bsub_fig  = 'bsub -q short -W 12:00 -R "rusage[mem=28000]" -o '
gatk_loc = '/groups/sorger/cchris/gatk/GenomeAnalysisTK.jar'
ref_fasta = '/groups/sorger/cchris/external_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

dir.create(outDir,showWarnings = F)
vcf_files <- grep('vcf$',list.files(vcfDir),value = T)

for (i in vcf_files){
  sample_id = gsub('.vcf','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out java -jar ',gatk_loc,' -T VariantFiltration -R ',ref_fasta,' -V ',vcfDir,i,' -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ',outDir,i,sep = ''))
}
