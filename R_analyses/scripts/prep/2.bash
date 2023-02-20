module load apps/R/4.0.3/gcc-8.2.0
Rscript 2.Integrate_normalise.R -d /fastdata/bop20pp/Avian_scRNAseq/R_analyses/outdata/filtered_seurat.RData -o /fastdata/bop20pp/Avian_scRNAseq/cellranger -t 1 \
	-c /home/bop20pp/software/MeioticDrive2022/R_analyses/data/cell_cycle_markers_complete.csv
