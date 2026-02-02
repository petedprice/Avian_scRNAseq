process seurat_count_matrix {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 84.GB * task.attempt }

    label 'seurat'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input:
    tuple val(species), file("${species}_integrated_seurat.RData"), file("${species}_motifs.tbl")

    output:
    tuple val(species), file("${species}_counts.csv"), file("${species}_mm_tfs.txt"), file("${species}_motifs.tbl")

    script:

    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/seurat_count_matrix.R \
        ${species} \
        ${species}_integrated_seurat.RData \
        ${species}_motifs.tbl \
        ${params.celltype_matrix}

    mv ${species}_seurat_counts.csv ${species}_counts.csv

    """
}
