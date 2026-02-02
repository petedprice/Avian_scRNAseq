process VariantFiltration {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    memory { 128.GB * task.attempt }


    input:
    tuple val(contig), file("${contig}_genotyped_allsites.vcf.gz")

    output:
    tuple val(contig), file("${contig}_genotyped_allsites_VF.vcf.gz")


    script:
    """
    #!/bin/bash
    gatk IndexFeatureFile -I ${contig}_genotyped_allsites.vcf.gz


    gatk --java-options "-Xmx4g" VariantFiltration \
	-V ${contig}_genotyped_allsites.vcf.gz \
	-O tmp.vcf.gz \
	--filter-expression "QUAL < 30.0" \
	--filter-name "LowQUAL" \
	--filter-expression "MQ < 40.0" \
	--filter-name "LowMQ" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "LowMQRankSum" \
	--filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "LowReadPosRankSum"

    mv tmp.vcf.gz ${contig}_genotyped_allsites_VF.vcf.gz


    """
}

