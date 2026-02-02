process SnpEff_anot {

    label 'snpf'
    errorStrategy 'retry'

    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(contig), file("SOMETHING.vcf.gz")

    output:
    tuple val(contig), file("${contig}_filt.recode.vcf")

    script:
    """
    #!/bin/bash

    snpEff.jar

    """
}

