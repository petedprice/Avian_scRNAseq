
process ctx {

    cpus { 60 * task.attempt }
    memory { 1999.GB * task.attempt }
    time = '48h'

    label 'SCENIC'

    publishDir 'regulons', mode: 'copy', overwrite: true, pattern: '*regulons.csv'

    input: 
    tuple val(species),  file("${species}_expr_mat.adjacencies.tsv"), file("${species}_motifs.tbl"), file("${species}_counts.csv"), file("combined")

    output:
    tuple val(species),  file("${species}_regulons.csv"), file("${species}_counts.csv")
    
    script:

    """
    #!/bin/bash
    cp combined/${species}_motif_db.regions_vs_motifs.rankings.feather .

    pyscenic ctx \
	${species}_expr_mat.adjacencies.tsv \
        ${species}_motif_db.regions_vs_motifs.rankings.feather \
        --annotations_fname ${species}_motifs.tbl \
        --output ${species}_regulons.csv \
        --num_workers $task.cpus \
        -t \
       --expression_mtx_fname ${species}_counts.csv 

    """
}


