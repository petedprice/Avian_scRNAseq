module load Nextflow/22.04.0
module load Anaconda3/2022.05

nextflow run /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/pop_gen/main.nf \
	--metadata /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/pop_gen/md.csv -resume \
	--fasta /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.fna \
	--gff /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.gff \
        --gtf /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.gtf \
	--pop_file /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/pop_gen/data/pop_file.txt \
	--protein /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.pro.faa \
	--cds /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.cds.fna \
	--assembly_report /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.AR.txt \
	--genbank /users/bi1pp/personal_git/Avian_scRNAseq/CL_analyses/nextflow/pop_gen/data/Anas_platyrhynchos_genbank.txt \
	--mrna /mnt/parscratch/users/bi1pp/Avian_DNA/refs/Anas_platyrhynchos.rna.fna



