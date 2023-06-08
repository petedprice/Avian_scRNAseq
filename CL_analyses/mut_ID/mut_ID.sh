#!/bin/bash

# Request one hour of wallclock time
#$ -l h_rt=4:0:0

# Request RAM
#$ -l rmem=1

# Select threads
#$ -pe smp 1


# Set the working directory
#$ -wd /fastdata/bop20pp/Avian_scRNAseq/mut_rate/wkdir

source /usr/local/extras/Genomics/.bashrc

source activate seurat

bash /home/bop20pp/software/Avian_scRNAseq/CL_analyses/mut_ID/mut_ID.bash $1 $2 $3
