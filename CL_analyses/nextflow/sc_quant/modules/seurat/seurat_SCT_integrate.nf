process seurat_SCT_integrate {

    cpus { 32 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 256.GB * task.attempt }    
    
    label 'seurat'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species),  val(stage), val(sex), file("${sample}_${species}_${sex}_${stage}_doublet_seurat.RData")

    output:
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_integrated_seurat.RData")

    script:
    """
    #!/bin/bash
    echo ${task.memory}

    Rscript ${projectDir}/Rscripts/seurat/SCT_integrate.R \
	. \
	${params.cellcycle_markers}

 	
    mv outdata/integrated_seurat.RData ${species}_${sex}_${stage}_integrated_seurat.RData



    """
}


