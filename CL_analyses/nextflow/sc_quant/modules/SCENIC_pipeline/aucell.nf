
process aucell {

    cpus { 60 * task.attempt }
    memory { 1999.GB * task.attempt }
    time = '48h'
    label 'SCENIC'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species),  file("${species}_regulons.csv"), file("${species}_counts.csv")

    output:
    tuple val(species),  file("${species}_regulons.csv")
    
    script:

    """
    #!/bin/bash
    pyscenic aucell \
        ${species}_counts.csv \
        ${species}_regulons.csv \
        -o ${species}_auc_mtx.csv \
        --num_workers $task.cpus \
        -t


    """
}


