process AddOrReplaceReadGroups {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    memory { 64.GB * task.attempt }

    input:
    tuple file("${sample}_sorted.bam"), file("${sample}_sorted.bam.bai"), val(sample)

    output:
    tuple file("${sample}_sorted.RG.bam"), val(sample)

    script:
    """
    #!/bin/bash

    RGID=@\$(samtools view ${sample}_sorted.bam | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
    RGLB="ILLUMINA_LB"
    RGSM=$sample
    RGPL="ILLUMINA"
    RGPU=\$RGID.\$RGPL

    gatk AddOrReplaceReadGroups \
	--I ${sample}_sorted.bam \
	--O ${sample}_sorted.RG.bam \
 	--RGID \$RGID \
 	--RGLB \$RGLB \
 	--RGPL \$RGPL \
 	--RGPU \$RGPU \
	--RGSM \$RGSM \
	--SORT_ORDER coordinate \
	--CREATE_INDEX TRUE

    
    """
}

