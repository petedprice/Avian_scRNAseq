
process seurat_count_matrix {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 84.GB * task.attempt }    
    
    label 'seurat'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_integrated_seurat.RData"), file("${species}_motifs.tbl")

    output:
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_counts.csv"), file("${species}_${sex}_${stage}_mm_tfs.txt"), file("${species}_motifs.tbl")

    script:

    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/seurat_count_matrix.R \
	${species} \
	${species}_${sex}_${stage}_integrated_seurat.RData \
	${species}_motifs.tbl

    mv ${species}_seurat_counts.csv ${species}_${sex}_${stage}_counts.csv
    mv ${species}_motifs.tbl ${species}_${sex}_${stage}_mm_tfs.txt

    """
}


