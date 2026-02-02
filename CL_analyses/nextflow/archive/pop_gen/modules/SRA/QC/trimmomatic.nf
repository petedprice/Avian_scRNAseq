process trimmomatic {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'java'


    input:
    tuple file('*fastq*'), val(SRA)

    output:
    tuple file("${SRA}_trimmed.fastq"), val(SRA)    
    
    script:
    """
    #!/bin/bash

    java -jar ${projectDir}/software/Trimmomatic-0.39/trimmomatic-0.39.jar \
	SE -phred33 \
	*fastq* \
	${SRA}_trimmed.fastq \
	SLIDINGWINDOW:4:20 MINLEN:25 \
	ILLUMINACLIP:${projectDir}/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:40:15

    """
}
