process grn {

    cpus { 60 * task.attempt }
    memory { 1999.GB * task.attempt }    
    time = '48h'    
    label 'SCENIC'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species), file("${species}_counts.csv"), file("${species}_mm_tfs.txt"), file("${species}_motifs.tbl")

    output:
    tuple val(species), file("${species}_expr_mat.adjacencies.tsv"), file("${species}_motifs.tbl"), file("${species}_counts.csv")
  
    script:

    """
    #!/bin/bash
    pyscenic grn \
        --num_workers $task.cpus \
        -o ${species}_expr_mat.adjacencies.tsv \
        -m grnboost2 -t \
        ${species}_counts.csv \
        ${species}_mm_tfs.txt


    """
}


