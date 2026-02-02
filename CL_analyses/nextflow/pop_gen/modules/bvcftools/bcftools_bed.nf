process bcftools_bed {

    label 'bcftools'
    errorStrategy 'retry'

    time = '8h'
    cpus { 1 * task.attempt }
    memory { 4.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:

    tuple val(contig), file("${contig}_filt.recode.vcf"), file(degen)

    output:
    tuple val(contig), file("fold_${contig}.vcf.gz"), env(bed_label)

    script:
    """
    #!/bin/bash
    bed_label=\$(basename "\$(readlink -f ${degen})" .bed)

    cat ${params.gff} | awk '{if (\$3 == "gene") print \$0;}' | grep ${contig} | cut -f1,4,5 > genes.bed

    cat $degen | grep $contig > tmp.degen.bed

    bcftools view ${contig}_filt.recode.vcf -Oz -o ${contig}_filt.recode.vcf.gz
    bcftools index ${contig}_filt.recode.vcf.gz


    bcftools view -R genes.bed ${contig}_filt.recode.vcf.gz -Oz -o tmp1.vcf.gz 
    bcftools index tmp1.vcf.gz
    bcftools view -R tmp.degen.bed tmp1.vcf.gz -Oz -o fold_${contig}.vcf.gz

    """
}

