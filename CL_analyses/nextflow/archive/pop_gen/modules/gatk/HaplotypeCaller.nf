process HaplotypeCaller {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    memory { 128.GB * task.attempt }

    //publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    tuple val(contig), file("${sample}_subset.bam"), val(sample)

    output:
    tuple val(contig), file("${sample}_${contig}.g.vcf.gz"), file("${sample}_${contig}.g.vcf.gz.tbi"), val(sample)

    script:
    """
    #!/bin/bash
    gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R ${params.ref} \
	-I ${sample}_subset.bam \
	-O ${sample}_${contig}.g.vcf.gz \
	-L $contig \
	-ERC GVCF

    

    """
}

