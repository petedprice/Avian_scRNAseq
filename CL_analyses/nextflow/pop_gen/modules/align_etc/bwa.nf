process bwa {
    cpus = 32
    memory = '64 GB'
    time = '12h'

    label 'bwa'

    input:
    tuple val(sample), file("${sample}_trim_reads")

    output:
    tuple val(sample), file("${sample}.sam")
    
    script:
    """
    #!/bin/bash
    ${projectDir}/software/bwa/bwa mem -M -t ${task.cpus} \
	${params.fasta} \
	${sample}_trim_reads/${sample}_R1_paired.fastq.gz \
	${sample}_trim_reads/${sample}_R2_paired.fastq.gz \
	> ${sample}.sam


    """
}
