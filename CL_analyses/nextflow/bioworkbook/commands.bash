module load apps/Nextflow/22.04.0/binary

nextflow run /home/bop20pp/software/Avian_scRNAseq/CL_analyses/nextflow/bioworkbook/main.nf \
	--fasta_dir /fastdata/bop20pp/Avian_scRNAseq/ref_files/ \
	--gtf_dir /fastdata/bop20pp/Avian_scRNAseq/ref_files/ \
	--cellranger /home/bop20pp/software/cellranger-7.0.1/cellranger \
	--metadata /home/bop20pp/software/Avian_scRNAseq/CL_analyses/nextflow/bioworkbook/metadata.csv \
	--read_dir /fastdata/bop20pp/Avian_scRNAseq/reads \
	--cellranger_data /fastdata/bop20pp/Avian_scRNAseq/cellranger \
	--run_cellranger FALSE \
	-resume \
	-with-dag flowchat.png \
	-with-trace
