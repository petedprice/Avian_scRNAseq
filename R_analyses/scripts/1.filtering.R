###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)


#filtering thresholds 
filt = "filtered"
ftr = 200

seurat_data <- Read10X(data.dir = 'data/filtered_feature_bc_matrix/')

seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                 min.features = ftr, ### FILTERING CELLS WITH LESS THAN 200 GENES
                                 project = 'test', 
                                 min.cells = 3) 

seurat_obj@meta.data$sample <- "test"


seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA) # Calculating the number of features per UMI 


#CALCULATING PROPORTION OF READS MAPPING TO MITOCHONDRIAL GENOME (x100 as originally a percentage)
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "mt")/ 100 
metadata <- seurat_obj@meta.data # Seperately creating metadata dataframe to save metrics etc without risking affecting seurat object
metadata$cells <- rownames(metadata) #add cell IDs to metadat
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
seurat_obj@meta.data <- metadata #Save the more complete metadat to the seurat object

#FILTERING DATA
filtered_seurat <- subset(x = seurat_obj, 
                           # (log10GenesPerUMI > 0.8) & # Can be dying cells or simple cell types such as blood cells
                            (mitoRatio < 0.05))
metadata_clean <- filtered_seurat@meta.data
filtered_seurat <- SCTransform(filtered_seurat, vars.to.regress = 
                                c("mitoRatio","nUMI"))
print("SelectIntegrationFeatures")

#Initial plot making for comparisons to after integration 

# Perform PCA
print("first unintegrated UMAPS")
seurat_SCT_normalised <- filtered_seurat
seurat_SCT_normalised <- RunPCA(object = seurat_SCT_normalised)
seurat_SCT_normalised <- RunUMAP(seurat_SCT_normalised, 
                                 dims = 1:30,
                                 reduction = "pca", 
                                 resolution = 0.2)

seurat_SCT_normalised <- FindNeighbors(seurat_SCT_normalised, dims = 1:30, verbose = FALSE)
seurat_SCT_normalised <- FindClusters(seurat_SCT_normalised, verbose = FALSE, resolution = 
                                        0.4)
DimPlot(seurat_SCT_normalised, label = T)
somatic = c(7)
germ = c(2,6,4,9)
FeaturePlot(seurat_SCT_normalised, features = 'DAZL')
germ_markers <- c("DAZL", "DDX4", "CLDN1", "DND1", "LOC101795050", "CARHSP1", "Vasa")
FeaturePlot(seurat_SCT_normalised, germ_markers)

somatic_markers <- c("POSTN", "LOC101798048", "SOX9", "DMRT1")

FeaturePlot(seurat_SCT_normalised, somatic_markers)
DimPlot(seurat_SCT_normalised, label = T)

seurat_SCT_normalised$cell_type <- "germ"
seurat_SCT_normalised$cell_type[seurat_SCT_normalised$SCT_snn_res.0.4 %in% somatic] <- 'somatic'
seurat_SCT_normalised$cell_type[seurat_SCT_normalised$SCT_snn_res.0.4 %in% germ] <- 'germ'

DimPlot(seurat_SCT_normalised, group.by ='cell_type')
