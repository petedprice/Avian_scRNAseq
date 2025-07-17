process seurat_doublet {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    //maxRetries 6
    memory { 84.GB * task.attempt }
    
    label 'seurat'

    //publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*doublet_seurat.RData'

    input: 
    tuple val(species), val(sample), val(stage), val(sex), file("filtered_seurat.RData"), file("${sample}_initial_QC_plots")

    output:
    tuple val(species),  file("*doublet_seurat.RData")

    script:
    """    
    #!/bin/bash

    Rscript ${projectDir}/Rscripts/seurat/doublet_finder.R \
	filtered_seurat.RData \
	. \
	${task.cpus} \
	${params.cellcycle_markers} \
	TRUE \
	${projectDir}
 	
    mv outdata/doublet_seurat.RData  ./${sample}_${species}_${sex}_${stage}_doublet_seurat.RData
    """

}


