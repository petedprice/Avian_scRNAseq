process fasterq_dump {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'sra_tools'

    input:
    val(SRA)

    output:
    tuple file('*fastq*'), val(SRA)
    
    
    script:
    """
    #!/bin/bash
    echo $SRA
    fasterq-dump $SRA 
    """
}
