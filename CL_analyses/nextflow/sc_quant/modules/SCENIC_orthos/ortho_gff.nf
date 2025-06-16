process ortho_gff {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'tidyverse'

    input:
    tuple file("all_genes.txt"), file("bird_species_orthos.txt"), file("hum_mouse_bird_species_orthos.txt"), file(gff1), file(gff2), file(gff3), file(gff4)

    output:
    file("OrthoDB_GFF_orthologs.csv")
    
    script:
    """
    #!/bin/bash
    Rscript ${projectDir}/Rscripts/SCENIC/ortho_gff.R ${launchDir}/libs


    """
}
