process seurat_soupX {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    //maxRetries 6
    memory { 84.GB * task.attempt }
    
    label 'seurat'


    input: 
    tuple val(species), val(sample), val(stage), val(sex), file("filtered_seurat.RData"), file("${sample}_crdata"), file("${sample}_initial_QC_plots")

    output:
    tuple val(species), val(sample), val(stage), val(sex), file("seurat_soupX.RData")

    script:
    """    
    #!/bin/bash


    Rscript ${projectDir}/Rscripts/seurat/soupX.R \
	filtered_seurat.RData \
	. \
	${task.cpus} \
	${params.cellcycle_markers} \
	${sample}_crdata/outs/raw_feature_bc_matrix



    """

}


