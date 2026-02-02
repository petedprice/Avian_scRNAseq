process seurat_doublet {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    //maxRetries 6
    memory { 84.GB * task.attempt }
    
    label 'seurat'

    input: 
    tuple val(species), val(sample), val(stage), val(sex), file("seurat_soupX.RData")

    output:
    tuple val(species),  val(stage), val(sex), file("${sample}_${species}_${sex}_${stage}_doublet_seurat.RData")

    script:
    """    
    #!/bin/bash

    Rscript ${projectDir}/Rscripts/seurat/doublet_finder_souped.R \
	seurat_soupX.RData \
	. \
	${task.cpus} \
	${params.cellcycle_markers} \
	TRUE \
	${projectDir}
 	
    mv outdata/doublet_seurat.RData  ./${sample}_${species}_${sex}_${stage}_doublet_seurat.RData
    """

}


