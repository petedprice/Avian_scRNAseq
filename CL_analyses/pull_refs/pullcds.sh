#!/bin/bash

#$ -l h_rt=4:0:0

#$ -l rmem=16G

#$ -pe smp 1

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/Avian_scRNAseq/wdir


rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ./
grep -E $1 assembly_summary_refseq.txt | cut -f 20 > ftp_links.txt
awk 'BEGIN{FS=OFS="/";filesuffix="cds_from_genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_cds_files.sh
source download_cds_files.sh
gunzip *.gz

#mv *cds* /fastdata/bop20pp/Avian_scRNAseq/ref_files

