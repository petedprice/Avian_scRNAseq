process R_mut_filtering {

    label 'tidyverse'

    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(species), val(sex), val(stage)

    output:
    tuple val(sample), file("${sample}_${contig}_snp_summarise.txt.gz")    
    script:
    """
    #!/bin/bash
    cat ${params.metadata} | grep $species | grep $sex | grep $stage > metadata_ss.csv
    Rscript ${projectDir}/Rscripts/seurat/1.filtering.R \
	-m ${projectDir}
	-d ${params.cellranger_data}
	-o .
	metadata_ss.csv
	$

    """
}
