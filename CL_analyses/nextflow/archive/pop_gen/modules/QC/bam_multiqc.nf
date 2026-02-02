process bam_multiqc {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'multiqc'

    publishDir 'multiqc', mode: 'copy', overwrite: true, pattern: 'bam_multi*'

    input:
    file("*_sorted_stats")

    output:
    file("bam_multi*")
    
    
    script:
    """
    #!/bin/bash
    multiqc .
    mv multiqc_report.html bam_multiqc_report.html
    mv multiqc_data bam_multiqc_data

    """
}
