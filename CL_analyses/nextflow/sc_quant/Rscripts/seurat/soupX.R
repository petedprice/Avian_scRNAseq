###LIBRARIES ----
library(Seurat)
library(Matrix)
library(ggplot2)

install.packages("SoupX", lib = '.')
library(SoupX, lib = '.')

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
load(args[1]) # path to filtered seurat RData
output_path <- args[2]
threads = as.numeric(args[3])-2
samples="all"
cellcycle <- read.csv(args[4], header = F)
memory=340
path_raw=args[5]

#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)


#### RUN COMMANDS ####

get_soup_groups <- function(sobj){
	sobj <-SCTransform(sobj, vars.to.regress =
                                c("nUMI","mitoRatio"))
	sobj <-RunPCA(sobj)
	sobj <-RunUMAP(sobj, dims = 1:40,reduction = "pca")
	sobj <-FindNeighbors(sobj, dims = 1:40,verbose = F)
	sobj <-FindClusters(sobj, verbose = F, resolution = 0.75)
	return(sobj@meta.data[['seurat_clusters']])
}
add_soup_groups <- function(sobj){
  sobj$soup_group <- get_soup_groups(sobj)
  return(sobj)
}

sobj <- add_soup_groups(filtered_seurat)


### IMPORT RAW MATRIX ###
raw.matrix <- Read10X(data.dir = path_raw)
rownames(raw.matrix) <- gsub("_", "-", rownames(raw.matrix))
filt.matrix <- GetAssayData(sobj, assay = "RNA", layer = "counts")
raw.matrix_subset <- raw.matrix[rownames(raw.matrix) %in% rownames(filt.matrix), ]  # Subset raw.matrix
sc <- SoupChannel(raw.matrix_subset, filt.matrix)
sc <- setClusters(sc, sobj$soup_group)
sc <- autoEstCont(sc, doPlot = F)
out <- adjustCounts(sc, roundToInt = TRUE)

#### SAVING COUNTS ####
sobj[["original.counts"]] <- CreateAssayObject(counts = filt.matrix)
sobj[["RNA"]] <- CreateAssayObject(counts = out)
sobj@meta.data$nCount_RNA <- colSums(out)

DefaultAssay(sobj) <- "RNA"

seurat_soupX <- sobj

save(seurat_soupX, file = "seurat_soupX.RData")
