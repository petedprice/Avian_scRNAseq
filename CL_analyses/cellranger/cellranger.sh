#!/bin/bash

#$ -l h_rt=8:0:0

#$ -l rmem=32G

#$ -pe smp 8

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/cellranger/wdir

name=$1
rn1='-'
id=${name#*$rn1}
id=${id%'/'}
echo $id

~/software/cellranger-7.0.1/cellranger count --id=$id \
	--fastqs=${2} \
	--transcriptome=${3}

