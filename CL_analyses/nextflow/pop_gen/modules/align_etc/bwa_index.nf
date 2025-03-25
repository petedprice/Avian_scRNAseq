process bwa_index {

    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }


    input:


    output:
    
    script:
    """
    #!/bin/bash

    ${projectDir}/software/bwa/bwa index ${params.fasta}


    """
}
