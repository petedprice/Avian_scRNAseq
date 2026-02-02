process CreateSequenceDictionary {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    memory { 128.GB * task.attempt }


    input:

    output:

    script:
    """
    #!/bin/bash
    gatk CreateSequenceDictionary -R ${params.fasta}

    samtools faidx ${params.fasta}

    """
}

