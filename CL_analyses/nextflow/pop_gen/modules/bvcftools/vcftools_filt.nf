process vcftools_filt {

    label 'vcftools'
    errorStrategy 'retry'

    time = '8h'
    cpus { 4 * task.attempt }
    memory { 16.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(contig), file("${contig}_genotyped_allsites_VF.vcf.gz")

    output:
    tuple val(contig), file("${contig}_filt.recode.vcf")

    script:
    """
    #!/bin/bash

    vcftools \
	--gzvcf ${contig}_genotyped_allsites_VF.vcf.gz \
	--out ${contig}_filt \
	--minGQ 30 \
        --minDP 10 \
        --recode \
        --remove-indels

    """
}

