process GenotypeGVCFs {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    memory { 128.GB * task.attempt }


    input:
    tuple val(contig), file("${contig}_database")

    output:
    tuple val(contig), file("${contig}_genotyped.vcf.gz")    

    script:
    """
    #!/bin/bash

    gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R ${params.ref} \
	-V gendb://${contig}_database \
	-L $contig \
	-all-sites \
	-O ${contig}_genotyped.vcf.gz

    """
}

