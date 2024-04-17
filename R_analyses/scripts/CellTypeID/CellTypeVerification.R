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

RDatas <- list.files("data/seurat_RData/", full.names = T)

sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"

markers <- read.csv("data/markers/markers.csv")

plot_PCA_and_markers <- function(x, type = "Feature"){
  
  load(x)
  print(x)
  #load(RDatas[1])
  DefaultAssay(seurat_marker) <- "RNA"
  
  si_tmp <- filter(sample_info, sample %in% unique(seurat_marker$sample))[1,]
  species <- si_tmp$species
  markers$marker <- markers[,species]
  save_name <- paste(si_tmp$species, si_tmp$sex, si_tmp$stage, sep = "_")
  plots <- list()
  markers$sex <- substring(markers$sex, 1, 1)
  mks_tmp <- filter(markers, sex %in% c("B", si_tmp$sex)) %>% 
    filter(marker %in% rownames(seurat_marker)) #%>% 
    #filter(celltype %in% unique(seurat_marker$sctype_labels))
  
  plots <- list()
  plots[[1]] <- DimPlot(seurat_marker, group.by = 'sctype_labels', label = T)  + NoLegend()
  plots[[2]] <- DimPlot(seurat_marker, group.by = 'seurat_clusters', label = T) + NoLegend()
  plots[[3]] <- DotPlot(seurat_marker, features = mks_tmp$marker, group.by = 'sctype_labels') +
    theme(axis.text.x=element_text(size=4),axis.text.y=element_text(size=7))

  for (i in 4:nrow(mks_tmp)){
    tmp_plots <- list()
    tmp_plots[[1]] <- FeaturePlot(seurat_marker, features = mks_tmp$marker[i]) + 
      ggtitle(mks_tmp$marker[i], mks_tmp$celltype[i])
    tmp_plots[[2]] <- VlnPlot(seurat_marker, features = mks_tmp$marker[i], group.by = 'sctype_labels')+ 
      ggtitle(mks_tmp$marker[i], mks_tmp$celltype[i])
    tmp_plots[[3]] <- VlnPlot(seurat_marker, features = mks_tmp$marker[i], group.by = 'seurat_clusters')+ 
      ggtitle(mks_tmp$marker[i], mks_tmp$celltype[i])
    
    plots <- c(plots, tmp_plots)
      
  }
  fin_plot <- DotPlot(seurat_marker, features = mks_tmp$marker, group.by = 'seurat_clusters')  +
    theme(axis.text.x=element_text(size=4),axis.text.y=element_text(size=7))
  plots <- c(plots, list(fin_plot))

  
  arranged <- ggarrange(plotlist = plots)
  ggsave(paste("plots/", type, "_", save_name, ".pdf", sep = ""), arranged, height = 49, width = 49)
}

#plot_PCA_and_markers(RDatas[1])
lapply(RDatas, plot_PCA_and_markers)
