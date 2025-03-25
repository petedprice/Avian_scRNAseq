library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(clustree)

sample_info <- read.csv("data/metadata_full.csv", header = F) 
colnames(sample_info) <- c('sample', "species", "ed", "sex", "folder", 'ref', "mito")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"
markers <- read.csv("data/markers/markers.csv")


RDatas <- list.files("data/seur_objs/integrated/", full.names = T, pattern = ".RData")
params <- data.frame(stage = NA, sex = NA, species = NA, res = NA, PCs = NA)
seur_objs <- list()

for (x in RDatas){
  print(x)
  load(x)
  DefaultAssay(seurat_integrated) <- "integrated"
  seurat_integrated <- RunPCA(seurat_integrated, verbose = F)
  
  #Decide resolution etc 
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

  seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:pcs, verbose = F)
  seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:pcs, verbose = F)
  
  print("MAPS RUN")
  ######### Determining resolution to use --------
  seurat_integrated <- FindNeighbors(object=seurat_integrated, dims=1:pcs, verbose = F)
  seurat_integrated <- FindClusters(object=seurat_integrated, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                             0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                                             2.5, 3), verbose = F)
  clusttree <- clustree(seurat_integrated, verbose = F)
  CT <- plot(clusttree)
  CT_noclusters <- clusttree$data %>% 
    group_by(integrated_snn_res.) %>% 
    summarise(no_clusters = n_distinct(cluster)) %>% 
    dplyr::rename(resolution = `integrated_snn_res.`) %>%
    ggplot(aes(x = resolution, y = no_clusters)) +
    geom_point()
  
  clustering_plots <- ggarrange(EP, CT, CT_noclusters, ncol = 1, nrow = 3)
  print("Clusters run, choose now")
  ggsave("plots/temp_clustering_plots.png", clustering_plots, height = 30, width = 30)
  system("open ./plots/temp_clustering_plots.png")
  res=NA
  while (!res %in% colnames(seurat_integrated@meta.data)){
    print("input your desired integration resolution")
    res=paste0('integrated_snn_res.', readline())
  }
  Idents(seurat_integrated) <- res
  
  system("rm plots/temp_clustering_plots.png")
  
  si_tmp <- filter(sample_info, sample %in% unique(seurat_integrated$sample))
  species <- si_tmp$species[1]
  sex = si_tmp$sex[1]
  stage = si_tmp$stage[1]
  
  seurat_integrated$sex <- sex
  seurat_integrated$stage <- stage
  seurat_integrated$species <- species
  
  plots <- list()
  mks_tmp <- filter(markers, sex %in% c("B", sex)) %>% 
    filter(marker %in% rownames(seurat_integrated))
  mks_tmp <- mks_tmp[,c("marker", "celltype", "sex", species)]
  plots <- list()
  plots[[1]] <- DimPlot(seurat_integrated, group.by = res, label = T) + NoLegend()
  plots[[2]] <- DotPlot(seurat_integrated, features = unique(mks_tmp$marker), group.by = res) + 
    coord_flip()
  arranged <- cowplot::plot_grid(plots[[1]], plots[[2]], nrow = 2, rel_heights = c(1,3))
  ggsave(paste("plots/numbered_clusters_", sex, "_", stage, "_", species, ".pdf", sep = ""), 
         arranged, height = 14, width = 5)
  
  params_temp <- data.frame(stage = stage, sex = sex, species = species, res = res, PCs = pcs)
  params <- rbind(params, params_temp)
  seur_objs[[paste(species, sex, stage, sep = "_")]] <- seurat_integrated
}

params <- params[!is.na(params$stage),]
#write.table(params, 'data/seur_objs/cluster_params.csv', sep = ',')

