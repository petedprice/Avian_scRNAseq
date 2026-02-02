process qualimap {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'qualimap'

    input:
    tuple file("${SRA}_sorted.bam"), file("${SRA}_sorted.bam.bai"), val(SRA)

    output:
    file("*_sorted_stats")
        
    script:
    """
    #!/bin/bash
    qualimap bamqc -bam ${SRA}_sorted.bam
    """
}
