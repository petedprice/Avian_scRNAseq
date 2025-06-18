module load Nextflow/22.04.0

nextflow run /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/sc_quant/main.nf \
	--cellranger /users/bi1pp/software/cellranger-9.0.0/cellranger \
	--metadata /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/sc_quant/metadata_apls.csv \
	--read_dir /mnt/parscratch/users/bi1pp/Avian_sc_quant/reads \
	--cellcycle_markers /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/sc_quant/data/chicken_cellcycle.csv \
	-resume
