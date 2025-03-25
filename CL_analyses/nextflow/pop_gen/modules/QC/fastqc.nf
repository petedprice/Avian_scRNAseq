process fastqc {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'fastqc'

    input:
    tuple file(reads), val(SRA)

    output:
    file("*_fastqc")
        
    script:
    """
    #!/bin/bash
    fastqc $reads 
    mkdir ${SRA}_fastqc
    mv  *fastqc.* ${SRA}_fastqc

    """
}
