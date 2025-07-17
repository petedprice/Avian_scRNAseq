
process ctx {

    cpus { 32 * task.attempt }
    memory { 128.GB * task.attempt }
    time = '8h'

    label 'SCENIC'

    publishDir 'regulons', mode: 'copy', overwrite: true, pattern: '*regulons.csv'

    input: 
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_expr_mat.adjacencies.tsv"), file("${species}_motifs.tbl"), file("${species}_${sex}_${stage}_counts.csv"), file("combined")

    output:
    tuple val(species),  val(stage), val(sex), file("${species}_${sex}_${stage}_regulons.csv"), file("${species}_${sex}_${stage}_counts.csv")
    
    script:

    """
    #!/bin/bash
    cp combined/${species}_motif_db.regions_vs_motifs.rankings.feather .

    pyscenic ctx \
	${species}_${sex}_${stage}_expr_mat.adjacencies.tsv \
        ${species}_motif_db.regions_vs_motifs.rankings.feather \
        --annotations_fname ${species}_motifs.tbl \
        --output ${species}_${sex}_${stage}_regulons.csv \
        --num_workers $task.cpus \
        -t \
       --expression_mtx_fname ${species}_${sex}_${stage}_counts.csv 

    """
}


