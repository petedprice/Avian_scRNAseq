process grn {

    cpus { 64 * task.attempt }
    memory { 2000.GB * task.attempt }    
    time = '48h'    
    label 'SCENIC'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_counts.csv"), file("${species}_${sex}_${stage}_mm_tfs.txt"), file("${species}_motifs.tbl")

    output:
    tuple val(species), val(stage), val(sex), file("${species}_${sex}_${stage}_expr_mat.adjacencies.tsv"), file("${species}_motifs.tbl"), file("${species}_${sex}_${stage}_counts.csv")
  
    script:

    """
    #!/bin/bash

    pyscenic grn \
        --num_workers $task.cpus \
        -o ${species}_${sex}_${stage}_expr_mat.adjacencies.tsv \
        -m grnboost2 -t \
        ${species}_${sex}_${stage}_counts.csv \
        ${species}_${sex}_${stage}_mm_tfs.txt
    echo done 

    """
}
