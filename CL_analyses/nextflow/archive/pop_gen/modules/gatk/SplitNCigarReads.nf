process SplitNCigarReads {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    publishDir 'dup_NCR_clean_bam', mode: 'copy', overwrite: true, pattern: '*dup_NCR.bam'

    input:
    tuple file("${sample}_sorted.RG.dups.bam"), val(sample)

    output:
    tuple file("${sample}_sorted.RG.dups.NCR.bam"), val(sample)

    script:
    """
    #!/bin/bash
    
    gatk SplitNCigarReads \
      -R ${params.ref} \
      -I ${sample}_sorted.RG.dups.bam \
      -O ${sample}_sorted.RG.dups.NCR.bam

    """
}

