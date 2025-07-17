process combine_partials {
    cpus = 4
    memory = '16 GB'
    time = '4h'

    label 'create_cistarget'

    input:
    tuple val(species), file(feathers)

    output:
    tuple val(species), file("combined")    

    script:
    """
    #!/bin/bash
    mkdir partials
    mkdir combined

    mv *feather* partials

    ${projectDir}/software/create_cisTarget_databases/combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py \
	-i partials \
	-o combined

     ${projectDir}/software/create_cisTarget_databases/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
	--db combined/${species}_motif_db.motifs_vs_regions.scores.feather

    """
}
