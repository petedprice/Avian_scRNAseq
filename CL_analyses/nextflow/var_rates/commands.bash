module load Nextflow/22.04.0
module load Anaconda3/2022.05

nextflow run /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/main.nf \
	--metadata /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/metadata_full.csv \
	--tree /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/data/tree.txt \
	--branch_trees /users/bop20pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/data/branch_trees \
	-resume



