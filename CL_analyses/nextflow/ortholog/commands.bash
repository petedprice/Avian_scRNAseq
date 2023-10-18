module load apps/Nextflow/22.04.0/binary

nextflow run /home/bop20pp/software/Avian_scRNAseq/CL_analyses/nextflow/ortholog/main.nf \
	--metadata /home/bop20pp/software/Avian_scRNAseq/CL_analyses/nextflow/ortholog/metadata_full.csv \
	--tree /home/bop20pp/software/Avian_scRNAseq/CL_analyses/nextflow/ortholog/data/tree.txt \
	-resume \
	-with-dag flowchat.png \
	-with-trace
