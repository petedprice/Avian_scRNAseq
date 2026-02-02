process bedtools_intersect {

    label 'bedtools'
    errorStrategy 'retry'

    time = '8h'
    cpus { 1 * task.attempt }
    memory { 4.GB * task.attempt }

    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:

    tuple val(contig), file("${contig}_filt.recode.vcf"), file(degen)

    output:
    tuple val(contig), file("fold_${contig}.vcf"), env(bed_label)

    script:
    """
    #!/bin/bash
    bed_label=\$(basename "\$(readlink -f ${degen})" .bed)

    cat ${params.gff} | awk '{if (\$3 == "gene") print \$0;}' | grep ${contig} | cut -f1,4,5 > genes.bed

    cat $degen | grep $contig > tmp.degen.bed

    bedtools intersect -a ${contig}_filt.recode.vcf -b genes.bed -wa > tmp1.vcf
    bedtools intersect -a	tmp1.vcf -b tmp.degen.bed -wa > fold_${contig}.vcf

    """
}

