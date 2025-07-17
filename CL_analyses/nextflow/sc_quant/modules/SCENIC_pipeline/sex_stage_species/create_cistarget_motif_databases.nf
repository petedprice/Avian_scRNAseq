process create_cistarget_motif_databases {
    cpus = 4
    memory = '16 GB'
    time = '4h'

    label 'create_cistarget'

    input:
    tuple file("v10nr_clust_public"), file("motifs.lst"), val(species), file("${species}_5000up_1000down_TSS.fasta"), val(part)    

    output:
    tuple val(species), file("${species}*.motifs_vs_regions.scores.feather")
    
    script:
    """
    #!/bin/bash

    ${projectDir}/software/create_cisTarget_databases/create_cistarget_motif_databases.py \
	-f ${species}_5000up_1000down_TSS.fasta \
	-M v10nr_clust_public/singletons \
	-m motifs.lst \
	-o ${species}_motif_db \
	-t $task.cpus \
	-p ${part} 250

    """
}
