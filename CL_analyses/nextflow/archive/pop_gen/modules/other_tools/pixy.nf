process pixy {

    label 'pixy'
    errorStrategy 'retry'


    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    input:
    tuple val(contig), file("${contig}_filt.recode.vcf")

    output:

    script:
    """
    #!/bin/bash

    echo $PATH

    cat ${params.gtf} | awk '{if (\$3 == "gene") print \$0;}' | grep ${contig} | cut -f1,4,5 > tmp.bed

    singularity run /mnt/parscratch/users/bi1pp/Avian_sc_quant/sing/pixy.sif bgzip ${contig}_filt.recode.vcf
    singularity run /mnt/parscratch/users/bi1pp/Avian_sc_quant/sing/pixy.sif tabix ${contig}_filt.recode.vcf.gz

    singularity run --bind .:/mnt/work \
	/mnt/parscratch/users/bi1pp/Avian_sc_quant/sing/pixy.sif \
	pixy --stats pi fst dxy \
	--vcf /mnt/work/${contig}_filt.recode.vcf.gz \
	--populations ${params.pop_file} \
	--n_cores ${task.cpus} \
	--bed_file /mnt/work/tmp.bed \
	--output_prefix /mnt/work/$contig


    """
}

