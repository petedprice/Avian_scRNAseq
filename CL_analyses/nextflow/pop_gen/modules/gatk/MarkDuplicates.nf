process MarkDuplicates {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    input:
    tuple file("${sample}_sorted.RG.bam"), val(sample)

    output:
    tuple file("${sample}_sorted.RG.dups.bam"), val(sample)

    script:
    """
    #!/bin/bash
    
    gatk MarkDuplicates \
	-I ${sample}_sorted.RG.bam \
	-O ${sample}_sorted.RG.dups.bam \
	-M ${sample}_dup_metrics.txt \
	--REMOVE_DUPLICATES

    """
}

