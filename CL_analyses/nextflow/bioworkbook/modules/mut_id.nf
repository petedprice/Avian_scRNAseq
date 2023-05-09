process mut_id {
    conda './envs/mut_id.yaml'

    input:
    tuple val(species), val(sample), val(contig), file("subset.bam"), file("${species}_cellranger_reference")

    output:

    script:
    """
    #!/bin/bash
    samtools index subset.bam
    
    ref=${species}_cellranger_reference/fasta/genome.fa
    
    #Filter bam for only barcoded reads
    samtools view -b -d CB -h -F 1024 -q 255 subset.bam > CB_subset.bam
    samtools index CB_subset.bam

    #Run BCFTOOLS
    bcftools mpileup --threads 1 -f \$ref CB_subset.bam | bcftools call --threads 1 -mv -Ob | bcftools view --types snps -i 'GT="het" &  INFO/DP >= 4 & (DP4[0]+DP4[1])>1 & (DP4[2]+DP4[3])>1' | gzip > snps.vcf.gz

    #RUN SAM2TSV
    java -jar ${projectDir}/software/jvarkit.ja sam2tsv -R \$ref --regions snps.vcf.gz -N -o sam2tsv.tsv.gz CB_subset.bam

    #GREP ALL SITES from bcftools in samttsv output to subset it
    zcat snps.vcf.gz | cut -f2 | egrep -v "^#"  | egrep -v POS > snps.txt
    zgrep -Ff snps.txt sam2tsv.tsv.gz  >  snps_sam2tsv.tsv
    cat snps_sam2tsv.tsv | cut -f1 | uniq > uniq_reads.txt

    #You now have read names supporting each variant in the file
    #You now need to associate the read names with the barcode.

    #Subset bam to just those reads
    samtools view CB_subset.bam | grep -Ff uniq_reads.txt > urCB_subset.sam

    # Extract barcodes from read names and pair with read name
    grep -o -E '.CB:Z.{0,19}' urCB_subset.sam | cut -f2 | cut -c 6- > barcodes.txt
    cut -f1 urCB_subset.sam > reads.txt
    paste reads.txt barcodes.txt  > reads_barcodes.txt

    #Merge reads_barcodes.txt with your snps_sam2tsv.tsv.gz file
    join reads_barcodes.txt snps_sam2tsv.tsv | gzip > ${sample}_${contig}_${species}_fin_snp_read.txt.gz

    """
}
