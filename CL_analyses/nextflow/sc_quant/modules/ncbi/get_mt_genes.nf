process get_mt_genes {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    input:
    tuple val(species), val(mt_contig), file("${species}.gtf"), file("${species}.fa")

    output:
    tuple val(species), file("${species}_mt_genes.txt")
    
    script:
    """
    #!/bin/bash

    cat ${species}.gtf | grep $mt_contig | cut -f9 | grep gene_id | sed -n 's/.*gene_id "\\([^"]*\\)".*/\\1/p' | uniq > ${species}_mt_genes.txt

    """
}
