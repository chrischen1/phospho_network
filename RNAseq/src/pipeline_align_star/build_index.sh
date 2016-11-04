module load seq/STAR/2.5.2a
out_dir = /groups/sorger/cchris/ref_genome
external_data = /groups/sorger/cchris/external_data
mkdir $out_dir
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $external_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $external_data/Homo_sapiens.GRCh38.85.gtf