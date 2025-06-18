
process ctx {

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 84.GB * task.attempt }    
    
    label 'SCENIC'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: '*RData'

    input: 
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_expr_mat.adjacencies.tsv"), file("${species}_motifs.tbl"), file("${species}_${sex}_${stage}_counts.csv"), ${species}_motif_db.regions_vs_motifs.rankings.feather

    output:
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_regulons.csv"), file("${species}_${sex}_${stage}_counts.csv")
    
    script:

    """
    #!/bin/bash
    pyscenic ctx \
	${species}_${sex}_${stage}_expr_mat.adjacencies.tsv
        ${species}_motif_db.regions_vs_motifs.rankings.feather \
        --annotations_fname ${species}_motifs.tbl \
        --output ${species}_${sex}_${stage}_regulons.csv \
        --num_workers ${nbr_threads} \
        -t \
       --expression_mtx_fname /data/${species}_seurat_counts.csv 

    """
}


