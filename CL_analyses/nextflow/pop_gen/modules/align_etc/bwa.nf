process bwa {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'bwa'

    input:
    tuple file("${SRA}_trimmed.fastq"), val(SRA)

    output:
    tuple file("${SRA}.sam"), val(SRA)

    
    script:
    """
    #!/bin/bash
    ${projectDir}/software/bwa/bwa mem -M -t 16 \
	${params.fasta} \
	${SRA}_trimmed.fastq \
	> ${SRA}.sam


    """
}
