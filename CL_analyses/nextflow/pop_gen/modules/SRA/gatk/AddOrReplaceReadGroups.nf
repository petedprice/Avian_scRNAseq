process AddOrReplaceReadGroups {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }

    input:
    tuple val(species), val(sample),  file("${sample}_crdata")

    output:
    tuple val(species), val(sample), file("${sample}_duplicates.bam")


    script:
    """
    #!/bin/bash

    RGID=$(zcat ${READS_1} | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4) # Extract RGID from read header
    RGLB="Illumina"
    RGSM=$sample
    RGPU=$RGID.$RGLB
    RGPL="Illumina"

    gatk AddOrReplaceReadGroups \
	--I 1_bams/${FILENAME}.sorted.bam \
	--O 1_bams/${FILENAME}.sorted.RG.bam \
 	--RGID $RGID \
 	--RGLB $RGLB \
 	--RGPL $RGPL \
 	--RGPU $RGPU \
	--RGSM $RGSM \
	--SORT_ORDER coordinate \
	--CREATE_INDEX TRUE

    
    """
}

