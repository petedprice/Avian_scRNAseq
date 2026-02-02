process trimmomatic {
    cpus = 16
    memory = '64 GB'
    time = '4h'

    label 'java'

    input:
    tuple val(sample), val(reads)

    output:
    tuple val(sample), file("${sample}_trim_reads")
    
    script:
    """
    #!/bin/bash

    java -jar ${projectDir}/software/Trimmomatic-0.39/trimmomatic-0.39.jar \
  	PE -phred33 \
  	${reads}/*_R1.fastq.gz \
  	${reads}/*_R2.fastq.gz \
  	${sample}_R1_paired.fastq.gz \
  	${sample}_R1_unpaired.fastq.gz \
  	${sample}_R2_paired.fastq.gz \
  	${sample}_R2_unpaired.fastq.gz \
  	ILLUMINACLIP:${projectDir}/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
  	SLIDINGWINDOW:4:20 \
  	MINLEN:50 \
	LEADING:3 \
	TRAILING:3 \
	-threads $task.cpus

    mkdir ${sample}_trim_reads
    mv ${sample}_*pair*gz ${sample}_trim_reads

    """
}
