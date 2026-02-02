process pixy {

    label 'pixy'
    errorStrategy 'retry'

    cpus { 1 * task.attempt }
    memory { 8.GB * task.attempt }

    publishDir 'pixy', mode: 'copy', overwrite: true, pattern: '*_pixy_out'

    input:
    tuple val(contig), file("fold_${contig}.vcf.gz"), val(bed_label)


    output:
    tuple val(contig), val(bed_label), path("${bed_label}_${contig}_pixy_out"), path(del)
    

    script:
    """
    #!/bin/bash

    cat ${params.gff} | awk '{if (\$3 == "gene") print \$0;}' | grep ${contig} | cut -f1,4,5 > genes.bed

    cp fold_${contig}.vcf.gz tmp.vcf.gz   
    #bgzip tmp.vcf
    tabix tmp.vcf.gz

    mkdir tmp_pixy_out
  
    pixy --stats pi tajima_d watterson_theta \
        --vcf tmp.vcf.gz \
        --populations ${params.pop_file} \
        --n_cores ${task.cpus} \
        --bed_file genes.bed \
        --output_folder ${bed_label}_${contig}_pixy_out


    cp -r ${bed_label}_${contig}_pixy_out del

    rm tmp.vcf.gz

    """
}

