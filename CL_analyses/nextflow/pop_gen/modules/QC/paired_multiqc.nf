process paired_multiqc {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'multiqc'

    tag {'paired_multiqc_' + '_' + sample }

    publishDir 'multiqc', mode: 'copy', overwrite: true, pattern: '*paired*'

    input:
    tuple val(sample), file('pre_fastqc'), file('post_fastqc')

    output:
    tuple val(sample), file("${sample}_paired_multiqc_data"), file("${sample}_paired_multiqc_report.html")
    
    script:
    """
    #!/bin/bash
    multiqc .
    mv multiqc_report.html ${sample}_paired_multiqc_report.html
    mv multiqc_data ${sample}_paired_multiqc_data

    """
}
