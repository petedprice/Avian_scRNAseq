
process aucell {

    cpus { 32 * task.attempt }
    memory { 128.GB * task.attempt }
    time = '8h'
    
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
        ${species}_${sex}_${stage}_counts.csv \
        ${species}_${sex}_${stage}_regulons.csv \
        -o ${species}_${sex}_${stage}_auc_mtx.csv \
        --num_workers $task.cpus \
        -t


    """
}


