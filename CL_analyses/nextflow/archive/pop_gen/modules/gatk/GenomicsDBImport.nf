process GenomicsDBImport {

    label 'gatk'
    cpus { 16 * task.attempt }
    memory { 128.GB * task.attempt }

    //publishDir 'raw_vcfs', mode: 'copy', overwrite: true, pattern: '*g.vcf.gz'

    input:
    tuple val(contig), file(vcfs), file(tbis), val(samples)

    output:
    tuple val(contig), file("${contig}_database")

    script:
    """
    #!/bin/bash
    ls -1 *.vcf.gz | sed 's/\\_${contig}\\.g\\.vcf\\.gz\$//' > tmp1
    ls -1 *vcf.gz > tmp2
    paste tmp1 tmp2 > cohort.sample_map

    gatk --java-options "-Xmx4g -Xms4g" \
	GenomicsDBImport \
	--genomicsdb-workspace-path ${contig}_database \
	-L ${contig} \
	--sample-name-map cohort.sample_map \
	--tmp-dir . \
	--reader-threads $task.cpus
    """
}

//    ls -1 *vcf.gz | sed 's/\.vcf\.gz\$//' > tmp1
