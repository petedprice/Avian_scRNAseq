process motif_modify {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'tidyverse'

    input:
    tuple file("v10nr_clust_public"), file("motifs.lst"), val(species), file("${species}_OrthoDB_GFF_orthologs.tsv")

    output:
    tuple val(species), file("${species}_motifs.tbl")    
    script:
    """
    #!/bin/bash


    Rscript ${projectDir}/Rscripts/SCENIC/modify_motif_collection.R ${species}_OrthoDB_GFF_orthologs.tsv $species

    """
}
