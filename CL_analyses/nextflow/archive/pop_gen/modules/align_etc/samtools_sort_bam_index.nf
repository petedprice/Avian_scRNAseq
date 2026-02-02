process samtools_sort_bam_index {
    cpus = 16
    memory = '64 GB'
    time = '4h'

    label 'samtools'

    input:
    val(sample)

    output:
    tuple file("${sample}_sorted.bam"), file("${sample}_sorted.bam.bai"), val(sample)
    
    script:
    """
    #!/bin/bash
    samtools view -b ${params.bams}/${sample}_crdata/outs/possorted_genome_bam.bam | samtools sort -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam

    """
}
