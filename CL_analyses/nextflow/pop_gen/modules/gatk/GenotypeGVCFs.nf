process GenotypeGVCFs {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 32 * task.attempt }
    memory { 512.GB * task.attempt }
    time = '8h'

    input:
    tuple val(contig), file("${contig}_database")

    output:
    tuple val(contig), file("${contig}_genotyped_allsites.vcf.gz")    

    script:
    """
    #!/bin/bash

    gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R ${params.fasta} \
	-V gendb://${contig}_database \
	-L $contig \
	-all-sites \
	-O ${contig}_genotyped_allsites.vcf.gz

    """
}

