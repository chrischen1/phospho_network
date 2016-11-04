fastq_path = '/groups/sorger/cchris/run_raw/'
out_path   = '/groups/sorger/cchris/run_qc/'
bsub_prefix = 'bsub -q short -W 4:00 -R "rusage[mem=18000]" '
trim_program = '/groups/sorger/cchris/Trimmomatic-0.36/trimmomatic-0.36.jar'


fastq_prefix <- unique(gsub('(\\w+)\\..+','\\1',grep('fastq',list.files(fastq_path),value = T)))
for (i in fastq_prefix){
  system(paste(bsub_prefix,' java -jar ',trim_program,' PE -phred33 ',
               fastq_path,i,'.R1.fastq ',fastq_path,i,'.R2.fastq ',
               ' -baseout ',out_path,i,'_qc.fastq ',
               'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36',sep = ''))
}