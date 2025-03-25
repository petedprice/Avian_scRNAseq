process samtools_sort_bam_index {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'samtools'

    input:
    tuple file("${SRA}.sam"), val(SRA)

    output:
    tuple file("${SRA}_sorted.bam"), file("${SRA}_sorted.bam.bai"), val(SRA)
    
    script:
    """
    #!/bin/bash
    samtools view -b ${SRA}.sam | samtools sort -o ${SRA}_sorted.bam
    samtools index ${SRA}_sorted.bam

    """
}
