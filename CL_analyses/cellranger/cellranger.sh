#!/bin/bash

#!/bin/bash

#$ -l h_rt=8:0:0

#$ -l rmem=16G

#$ -pe smp 4

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/Avian_scRNAseq/wdir

name=$1
rn1='-'
id=${name#*$rn1}
id=${id%'/'}
echo $id

~/software/cellranger-7.0.1/cellranger count --id=$id \
	--fastqs=${1}/ \
	--transcriptome=${2}

cp $id /fastdata/bop20pp/Avian_scRNAseq/cellranger
