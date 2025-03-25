#!/usr/bin/env Rscript
###LIBRARIES ----
library(Seurat)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggplot2)
library(future)
library(parallel)

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
cellcycle <- read.csv(args[2], header = F)
colnames(cellcycle) <- c("gene", "phase")
output_path <- args[1]
#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)


seurat_objs <- list.files(pattern = ".RData")
split_seurat <-list()
for (sobj in seurat_objs){
	load(sobj)
	print("SCT transform")
        doublet_seurat <- NormalizeData(doublet_seurat)
	doublet_seurat <- CellCycleScoring(doublet_seurat, 
                                    g2m.features = cellcycle$gene[cellcycle$phase == "G2/M"], 
                                    s.features = cellcycle$gene[cellcycle$phase == "S"])

	doublet_seurat <- SCTransform(doublet_seurat, vars.to.regress =
                                c("nUMI","mitoRatio"))
	split_seurat[[doublet_seurat$sample[1]]] <-doublet_seurat

}



#### #prep data for integration ---- 
# Identify variable features for integrating
print("SelectIntegrationFeatures")
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

#Prepossessing step necessary if SCT transformed
print("PrepSCTIntegration")
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = features)

#Find anchors that link datasets
print("FindIntegrationAnchors")
anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                  anchor.features = features, 
                                  normalization.method = "SCT")

#Initial plot making for comparisons to after integration 

##### first Perform PCA -----
print("first unintegrated UMAPS")
seurat_SCT_normalised <- Reduce(merge, split_seurat)
seurat_SCT_normalised <- RunPCA(object = seurat_SCT_normalised, features = features)
seurat_SCT_normalised <- RunUMAP(seurat_SCT_normalised, 
                                 dims = 1:30,
                                 reduction = "pca")
gc()
seurat_SCT_normalised <- FindNeighbors(seurat_SCT_normalised, dims = 1:30, verbose = FALSE)
seurat_SCT_normalised <- FindClusters(seurat_SCT_normalised, verbose = FALSE)

save(seurat_SCT_normalised, file = paste(outdatapath, "/seurat_SCT_normalised.RData", sep = ""))

fs_PCA1 <- DimPlot(seurat_SCT_normalised,
                   split.by = "sample")
ggsave(filename = paste(plotpath, "fs_sample_PCA.pdf", sep = ""))
fs_PCA2 <- DimPlot(seurat_SCT_normalised,
                   split.by = "Phase", group.by = "Phase")
ggsave(filename = paste(plotpath, "fs_Phase_PCA.pdf", sep = ""))


gc()
#Integrate data 
print("integrating")
seurat_integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT")

print("running PCA")
seurat_integrated <- RunPCA(object = seurat_integrated)

EP <- ElbowPlot(seurat_integrated, ndims = 40)
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
  
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
  
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
# last point where change of % of variation is more than 0.1%.
pcs <- min(co1, co2)
print("PCS decided")

seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:pcs,
                             reduction = "pca")

print("running UMAP")
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:pcs,
                             reduction = "pca")


seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:pcs)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                                 0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                 2.5, 3))
#Save data
print("saving data")
save(seurat_integrated, file = paste(outdatapath, "/integrated_seurat.RData", sep = ""))
save(split_seurat, anchors, features, file = paste(outdatapath, "/integrated_seurat_extras.RData", sep = ""))

