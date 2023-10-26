module load Nextflow/22.04.0

nextflow run /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/ortholog/main.nf \
	--metadata /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/ortholog/metadata_full.csv \
	--tree /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/ortholog/data/tree.txt \
        --paml_tree /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/ortholog/data/paml_tree.txt \
	-resume \
	-with-dag flowchat.png \
	-with-trace




