
process aucell {

    cpus { 32 * task.attempt }
    memory { 256.GB * task.attempt }
    time = '48h'


    label 'SCENIC'

    publishDir 'aucells', mode: 'copy', overwrite: true, pattern: '*mtx.csv'

    input: 
    tuple val(species), val(stage), val(sex), file("${species}_${sex}_${stage}_regulons.csv"), file("${species}_${sex}_${stage}_counts.csv")
//    tuple val(species),  file("${species}_regulons.csv"), file("${species}_counts.csv")

    output:
    tuple val(species), val(stage), val(sex), file("${species}_${sex}_${stage}_auc_mtx.csv")
    
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


