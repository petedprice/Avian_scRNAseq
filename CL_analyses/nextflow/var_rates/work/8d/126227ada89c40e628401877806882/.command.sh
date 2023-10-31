#!/bin/bash
echo Numida_meleagris
echo GCF_002078875.1
rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ./
grep -E GCF_002078875.1 assembly_summary_refseq.txt | cut -f 20 > ftp_links.txt
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_gff_files.sh
awk 'BEGIN{FS=OFS="/";filesuffix="cds_from_genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_cds_files.sh
awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz "}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_pro_files.sh

source download_gff_files.sh
source download_cds_files.sh
source download_pro_files.sh

gunzip *gz
mv *genomic.gff Numida_meleagris.gff
mv *cds_from_genomic.fna Numida_meleagris.cds.fna
mv *protein.faa Numida_meleagris.protein.faa
