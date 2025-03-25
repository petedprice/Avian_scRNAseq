process multiqc {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'multiqc'

    publishDir 'multiqc', mode: 'copy', overwrite: true, pattern: 'multiqc*'

    input:
    file("*_fastqc")

    output:
    file("multiqc*")   
    
    script:
    """
    #!/bin/bash
    multiqc .

    """
}
