process seurat_mutation {

    label 'seurat'

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }

    publishDir 'seurat_RData', mode: 'copy', overwrite: true, pattern: '*_seurat_mutant.RData'

    input: 
    tuple val(samples), val(sex), val(stage), val(species), file(mutation_data), file("marker_seurat.RData")


    output:
    tuple val(samples), val(sex), val(stage), val(species), file(mutation_data), file("${species}_${sex}_${stage}_seurat_mutant.RData")

    script:
    """
    #!/bin/bash

    cat ${params.metadata} | grep $species | grep ,$sex, | grep $stage > metadata_ss.csv

    echo $samples | sed 's/[],[]//g' > samples.txt
    Rscript ${projectDir}/Rscripts/seurat/4.seurat_mutation.R \
	marker_seurat.RData \
	. \
	samples.txt \
	metadata_ss.csv
    
    mv outdata/seurat_mutant.RData ${species}_${sex}_${stage}_seurat_mutant.RData 


    """
}


