process ortho_gff {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'tidyverse'

    input:
    tuple val(species), file("all_genes.txt"), file("hum_mouse_${species}_chicken_orthos.txt"), file("${species}_chicken_orthos.txt"), file("Chicken.gtf"), file("${species}.gtf")

    output:
    tuple val(species), file("${species}_OrthoDB_GFF_orthologs.tsv")
    
    script:
    """
    #!/bin/bash




    Rscript ${projectDir}/Rscripts/SCENIC/ortho_gff.R ${launchDir}/libs ${params.metadata} ${species}_chicken_orthos.txt hum_mouse_${species}_chicken_orthos.txt $species
    mv OrthoDB_GFF_orthologs.tsv ${species}_OrthoDB_GFF_orthologs.tsv



    """
}
