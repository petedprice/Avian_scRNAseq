process qualimap {
    cpus = 16
    memory = '64G'
    time = '8h'

    label 'qualimap'

    input:
    tuple file("${sample}_sorted.bam"), file("${sample}_sorted.bam.bai"), val(sample)

    output:
    file("*_sorted_stats")
        
    script:
    """
    #!/bin/bash
    qualimap bamqc -nt ${task.cpus} --java-mem-size=${task.memory} -bam ${sample}_sorted.bam
    """
}
