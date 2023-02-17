#!/bin/bash

#$ -l h_rt=4:0:0

#$ -l rmem=16G

#$ -pe smp 4

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/Avian_scRNAseq/wdir

~/software/cellranger-7.0.1/cellranger mkref \
	--genome=${3} \
	--fasta=${1} \
	--genes=${2}

mv $3 /fastdata/bop20pp/Avian_scRNAseq/ref_files/cellranger

