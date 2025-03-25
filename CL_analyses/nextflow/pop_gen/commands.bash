module load Nextflow/22.04.0
module load Anaconda3/2022.05

nextflow run /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/pop_gen/main.nf \
	--metadata /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/pop_gen/metadata_full.csv \
	--fasta /mnt/parscratch/users/bi1pp/Avian_pop_gen/refs/Mallard.fa \
	--gtf /mnt/parscratch/users/bi1pp/Avian_pop_gen/refs/Mallard.gtf \
	-resume



