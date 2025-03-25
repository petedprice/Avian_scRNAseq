process trimmed_multiqc {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'multiqc'

    publishDir 'multiqc', mode: 'copy', overwrite: true, pattern: 'trimmed*'

    input:
    file("*_fastqc")

    output:
    file("trimmed*")
    
    
    script:
    """
    #!/bin/bash
    multiqc .
    mv multiqc_report.html trimmed_multiqc_report.html
    mv multiqc_data trimmed_multiqc_data

    """
}
