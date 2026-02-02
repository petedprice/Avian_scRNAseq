process vcftools_filt {

    label 'vcftools'
    errorStrategy 'retry'

    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(contig), file("${contig}_genotyped.vcf.gz")

    output:
    tuple val(contig), file("${contig}_filt.recode.vcf")

    script:
    """
    #!/bin/bash

    vcftools \
	--gzvcf ${contig}_genotyped.vcf.gz \
	--out ${contig}_filt \
	--minGQ 30 \
	--minDP 10 \
	--minQ 500 \
	--max-missing 0.5 \
	--max-alleles 2 \
	--recode \
	--remove-indels \
	--recode-INFO-all



    """
}

