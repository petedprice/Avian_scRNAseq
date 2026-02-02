process split_bam_contig {
    
    cpus = 8
    memory = '32 GB'
    time = '4h' 


    label 'gatk' 
    
    input:
    tuple file("${sample}_sorted.RG.dups.bam"), val(sample), file("contig.txt")

    output:
    tuple env(contig), file("${sample}_subset.bam"), val(sample)

    script: 
    """    
    samtools index ${sample}_sorted.RG.dups.bam

    contig=\$(cat contig.txt)
    samtools view -bh ${sample}_sorted.RG.dups.bam \$contig > ${sample}_subset.bam
    samtools index ${sample}_subset.bam

    """
}
