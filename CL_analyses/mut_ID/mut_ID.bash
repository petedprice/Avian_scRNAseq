ref=$1
sample=$2

#Filter bam for only barcoded reads
samtools view -b -d CB ${sample}.bam > CB_sample.bam 
samtools index CB_sample.bam 

#Run BCFTOOLS 
bcftools mpileup --threads 1 -f $ref CB_sample.bam | bcftools call --threads 1 -mv -Ob | bcftools view -i 'GT="het"' | bcftools view -i "INFO/DP >= 3" | gzip > snps.vcf.gz

#RUN SAM2TSV 
java -jar ~/software/jvarkit/dist/jvarkit.jar sam2tsv -R $ref --regions snps.vcf.gz -N -o sam2tsv.tsv.gz CB_sample.bam

#GREP ALL SITES from bcftools in samttsv output to subset it 
zcat snps.vcf.gz | cut -f2 | egrep -v "^#"  | egrep -v POS > snps.txt
zgrep -Ff snps.txt sam2tsv.tsv.gz  >  snps_sam2tsv.tsv
cat snps_sam2tsv.tsv | cut -f1 | uniq > uniq_reads.txt

#You now have read names supporting each variant in the file 
#You now need to associate the read names with the barcode. 

# Extract barcodes from read names and pair with read name
samtools view CB_sample.bam | grep -Ff uniq_reads.txt | grep -o -E '.CB.{0,22}' | cut -d: -f3 > barcodes.txt
paste uniq_reads.txt barcodes.txt  > reads_barcodes.txt

#Merge reads_barcodes.txt with your snps_sam2tsv.tsv.gz file 
join reads_barcodes.txt snps_sam2tsv.tsv | gzip > fin_snp_read.txt.gz

