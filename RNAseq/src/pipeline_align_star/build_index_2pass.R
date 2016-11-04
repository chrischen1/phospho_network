gff_file = '/groups/sorger/cchris/external_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
gtf_file = '/groups/sorger/cchris/external_data/Homo_sapiens.GRCh38.85.gtf'
sjdb_path = '/groups/sorger/cchris/bam_1pass/'
outDir    = '/groups/sorger/cchris/ref_genome_2pass/'
bsub_fig  = 'bsub -q long -W 72:00 -R "rusage[mem=58000]" -o ref_index_2pass.out'

system('module load seq/STAR/2.5.2a')
dir.create(outDir,showWarnings = F)
sjdb_files <- paste0(paste(sjdb_path,grep('SJ',list.files(sjdb_path),value = T),sep = ''),collapse = ' ')

system(paste(bsub_fig,' STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ',outDir,' --genomeFastaFiles ',gff_file,' --sjdbGTFfile ',gtf_file,' --sjdbOverhang 75  --sjdbFileChrStartEnd ',sjdb_files,sep = ''))