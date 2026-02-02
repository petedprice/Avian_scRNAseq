module load Nextflow/22.04.0
module load Anaconda3/2022.05

nextflow run /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/main.nf \
	--metadata /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/WF_metadata.csv \
	--tree /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/data/WF_tree.txt \
	--branch_trees /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/data/WF_branch_trees \
	-resume



