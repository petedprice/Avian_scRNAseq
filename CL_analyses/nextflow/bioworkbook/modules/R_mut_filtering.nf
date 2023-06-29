process R_mut_filtering {

    label 'tidyverse'

    cpus { 18 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 250.GB * task.attempt }


    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(species), val(sample), val(contig), file("${species}_${sample}_${contig}_fin_snp_read.txt.gz")

    output:
    tuple val(sample), file("${sample}_${contig}_snp_summarise.txt.gz")    
    script:
    """
    #!/bin/bash
    Rscript ${projectDir}/Rscripts/mut_filtering.R \
	${species}_${sample}_${contig}_fin_snp_read.txt.gz \
	$sample \
	$contig 
    echo printing

    """
}
