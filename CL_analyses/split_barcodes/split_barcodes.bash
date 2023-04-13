export BAM_FILE='/fastdata/bop20pp/Avian_scRNAseq/cellranger/apls_12/outs/possorted_genome_bam.bam'
samtools view -H $BAM_FILE > SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f filter.txt > filtered_SAM_body

