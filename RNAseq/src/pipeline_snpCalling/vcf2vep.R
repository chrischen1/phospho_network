vcfDir  = '/groups/sorger/cchris/vcf_files_filtered/'
outDir    = '/groups/sorger/cchris/vep_files/'
bsub_fig  = 'bsub -q short -W 2:00 -o '
vep       = 'perl ~/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl'

dir.create(outDir,showWarnings = F)
vcf_files <- grep('.vcf$',grep('vep',list.files(vcfDir),value = T,invert = T),value = T)

for (i in vcf_files){
  system('module load seq/vep/83')
  sample_id = gsub('.vcf','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out ',vep,' -i ',vcfDir,i,' -offline -symbol -o ',outDir,sample_id,'.vep.txt',sep = ''))
}