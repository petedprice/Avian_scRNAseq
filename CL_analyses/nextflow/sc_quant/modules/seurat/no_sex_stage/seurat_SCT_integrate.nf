process seurat_SCT_integrate {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 256.GB * task.attempt }    
    
    label 'seurat'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species), file("doublet_seurat.RData")

    output:
    tuple val(species),  file("${species}_integrated_seurat.RData")

    script:
    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/SCT_integrate.R \
	. \
	${params.cellcycle_markers}

 	
    mv outdata/integrated_seurat.RData ${species}_integrated_seurat.RData


    """
}


