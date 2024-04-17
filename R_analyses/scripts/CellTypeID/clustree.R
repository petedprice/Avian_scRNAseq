library(SCINA)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(glmmTMB)
library(lme4)
library(clustree) #####

load("data/seurat_RData/duck_F_ED17_marker_seurat.RData")

sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"

markers <- read.csv("data/markers/markers.csv")


##################### CLUSTREE -------------
clustree(seurat_marker, prefix = "integrated_snn_res.")

########### JACKSTRAW -------------
seurat_marker <- JackStraw(seurat_marker, num.replicate = 100)
JackStrawPlot(seurat_marker, dims = 1:15)


######## HBC TRAINING -----
ElbowPlot(seurat_marker, ndims = 40)
colums <- colnames(seurat_marker@meta.data)
clusters <- seurat_marker@meta.data[,colums[startsWith(colums, "inte")]] %>% 
  as.data.frame() %>% 
  apply(., 2, function(x)(return(max(as.numeric(x)))))

resolution <- which.min(abs(clusters - 20)) %>% names()
seurat_marker[['new_clusters']] <- seurat_marker[[resolution]]

table(seurat_marker[['new_clusters']])
