
process aucell {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 84.GB * task.attempt }    
    
    label 'SCENIC'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_regulons.csv"), file("${species}_${sex}_${stage}_counts.csv")

    output:
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_regulons.csv")
    
    script:

    """
    #!/bin/bash
    pyscenic aucell \
        ${species}_${sex}_${stage}_seurat_counts.csv \
        ${species}_${sex}_${stage}_regulons.csv \
        -o ${species}_${sex}_${stage}_auc_mtx.csv \
        --num_workers 1
        -t


    """
}


