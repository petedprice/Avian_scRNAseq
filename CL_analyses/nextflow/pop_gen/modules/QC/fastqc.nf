process fastqc {
    cpus = 8
    memory = '32 GB'
    time = '4h'

    label 'fastqc'

    input:
    tuple val(sample), val(reads)

    output:
    tuple val(sample), file("*_fastqc")
        
    script:
    """
    #!/bin/bash
    fastqc ${reads}/*fastq.gz -o .
    mkdir ${sample}_fastqc
    mv  *fastqc.* ${sample}_fastqc

    """
}
