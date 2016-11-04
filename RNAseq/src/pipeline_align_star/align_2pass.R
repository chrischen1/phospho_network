genomeDir = '/groups/sorger/cchris/ref_genome_2pass'
fastqDir  = '/groups/sorger/cchris/run_qc/'
outDir    = '/groups/sorger/cchris/bam_2pass/'
bsub_fig  = 'bsub -q short -W 4:00 -R "rusage[mem=38000]" -o '

system('module load seq/STAR/2.5.2a')
dir.create(outDir,showWarnings = F)
fastq_prefix <- unique(gsub('(\\w+)_1P\\.fastq','\\1',grep('1P.fastq',list.files(fastqDir),value = T)))

for (i in fastq_prefix){
  file1 <- paste(fastqDir,i,'_1P.fastq',sep = '')
  file2 <- paste(fastqDir,i,'_2P.fastq',sep = '')
  system(paste(bsub_fig,outDir,i,'.out ','STAR --runThreadN 2 --genomeDir ',genomeDir,' --readFilesIn ',file1,' ',file2,' --outFileNamePrefix ',outDir,i,sep = ''))
}