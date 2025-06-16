process modify_motif_collection {
    cpus = 1
    memory = '4 GB'
    time = '4h'
 
    label 'tidyverse'

    input:
    tuple file("motifs.lst"), file("v10nr_clust_public")

    output:

    
    script:
    """
    #!/bin/bash

    Rscript ${projectDir}/Rscripts/SCENIC/modify_motif_collection.R v10nr_clust_publicv

    """
}
