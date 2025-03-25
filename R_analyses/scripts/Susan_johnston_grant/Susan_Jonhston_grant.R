library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(lme4)
library(RColorBrewer)

sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"
markers <- read.csv("data/markers/markers.csv")
clusters_df <- readxl::read_excel("data/Susan_Johnston_Grant_Figures/SusanJohnston_Grant_CelltypeIDs.xlsx", sheet = 1)
load_func <- function(x){
  load(x)
  return(seurat_marker)
}

RData1 <- list.files("data/seurat_RData/", full.names = T, pattern = c("ED24"))
RData2 <- list.files("data/seurat_RData/", full.names = T, pattern = c("ED22"))
RData3 <- list.files("data/seurat_RData/", full.names = T, pattern = c("ED25"))
RDatas <- c(RData1, RData2, RData3)

seur_objs <- lapply(RDatas, load_func)
length(seur_objs)
names(seur_objs) <- str_split(RDatas, "/", simplify = T)[,4] %>% 
  gsub("_marker_seurat.RData", "", .)


plot_PCA_and_markers <- function(x, type = "Feature"){
  seurat_marker <- x
  DefaultAssay(seurat_marker) <- "RNA"

  si_tmp <- filter(sample_info, sample %in% unique(seurat_marker$sample))
  species <- si_tmp$species
  plots <- list()
  markers$sex <- substring(markers$sex, 1, 1)
  mks_tmp <- filter(markers, substring(sex, 1, 1) %in% c("B", si_tmp$sex)) %>% 
    filter(marker %in% rownames(seurat_marker))
  mks_tmp <- mks_tmp[,c("marker", "celltype", "sex", species[1])]
  plots <- list()
  plots[[1]] <- DimPlot(seurat_marker, group.by = 'integrated_snn_res.0.1', label = T) + NoLegend()
  plots[[2]] <- DotPlot(seurat_marker, features = unique(mks_tmp$marker), group.by = 'integrated_snn_res.0.1') + 
    coord_flip()
  arranged <- cowplot::plot_grid(plots[[1]], plots[[2]], nrow = 2, rel_heights = c(1,3))
  ggsave(paste("plots/numbered_clusters", si_tmp$sex[1], "_", si_tmp$stage[1], "_", species[1], ".pdf", sep = ""), arranged, height = 14, width = 5)
  
}


celltypes <- clusters_df[,-1] %>% unlist() %>% unique()
celltypes <- celltypes[!celltypes %in% c("NA", "NO_CLUSTER")]
myColors <- brewer.pal(n = 12, "Paired")[-2]
names(myColors) <- levels(celltypes) 

lapply(seur_objs, plot_PCA_and_markers)
x <- seur_objs[[1]]


plot_PCA_and_markers2 <- function(x, type = "Feature"){
  seurat_marker <- x
  DefaultAssay(seurat_marker) <- "RNA"
  
  si_tmp <- filter(sample_info, sample %in% unique(seurat_marker$sample))
  sex = si_tmp$sex
  species <- si_tmp$species
  
  cl_tmp <- clusters_df[,c("Cluster", paste(species[1], sex[1], sep = " "))]
  clusters <- seurat_marker$integrated_snn_res.0.1
  new_clusters <- cl_tmp[,2] %>% unlist()
  names(new_clusters) <- cl_tmp$Cluster %>% unlist()
  
  new_clusters2 <- new_clusters[clusters]
  new_clusters2
  seurat_marker@meta.data$susan_clusters <- new_clusters2
  
  seurat_marker <- subset(seurat_marker, susan_clusters != "NA")
  seurat_marker@meta.data <- 
    seurat_marker@meta.data %>% 
    mutate(susan_clusters = factor(susan_clusters, levels = celltypes))
         
  plots <- list()
  markers$sex <- substring(markers$sex, 1, 1)
  mks_tmp <- filter(markers, substring(sex, 1, 1) %in% c("B", si_tmp$sex)) %>% 
    filter(marker %in% rownames(seurat_marker))
  mks_tmp <- mks_tmp[,c("marker", "celltype", "sex", species[1])]
  plots <- list()
  sex = sex[1]
  stage = si_tmp$ed[1]
  if (sex == "F"){
    sex = "female"
  } else { 
    sex = "male"}
  plots[[1]] <- DimPlot(seurat_marker, group.by = 'susan_clusters', label = F) + NoLegend() + 
    ggtitle(paste(sex, species, ", stage:" , stage, sep = " ")) + 
    labs(x = "UMAP-1", y = "UMAP-2") + 
    theme_classic() + 
    scale_colour_manual(values=myColors)
  
  plots[[2]] <- DotPlot(seurat_marker, features = unique(mks_tmp$marker), group.by = 'susan_clusters') + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
          axis.ticks = element_line(color = "black")) + 
    labs(x = "", y = "")
  cell_numbers <- table(seurat_marker$susan_clusters) %>% as.data.frame()
  cell_numbers$sex = sex
  cell_numbers$species = species[1]
  colnames(cell_numbers)[1] <- "Celltype"
  plots[[3]] <- cell_numbers
  arranged <- cowplot::plot_grid(plots[[1]], plots[[2]], nrow = 2, rel_heights = c(1.8,3))
  return(plots)
  ggsave(paste("plots/celltype_clusters_", si_tmp$sex[1], "_", si_tmp$stage[1], "_", species[1], ".pdf", sep = ""), arranged, height = 14, width = 7)
}

plots <- lapply(seur_objs, plot_PCA_and_markers2)

A <- plots$duck_F_ED24[[1]] + theme(legend.position = c(0.2, 0.2))
B <- plots$pheasant_F_ED22[[1]] + guides(colour = 'none')
C <- plots$guineafowl_F_ED25[[1]] + guides(colour = 'none')
D <- plots$duck_F_ED24[[2]]

figure <- cowplot::plot_grid(A, B, C, D, labels = c("A", "B", "C", "D"))
ggsave("data/Susan_Johnston_Grant_Figures/Females_figure.pdf", height = 15, width = 15)
system("open data//Susan_Johnston_Grant_Figures/Females_figure.pdf")


A <- plots$duck_M_ED24[[1]] + theme(legend.position = c(0.8, 0.2))
B <- plots$pheasant_M_ED22[[1]] + guides(colour = 'none')
C <- plots$guineafowl_M_ED25[[1]] + guides(colour = 'none')
D <- plots$duck_M_ED24[[2]]

figure <- cowplot::plot_grid(A, B, C, D, labels = c("A", "B", "C", "D"))
ggsave("data/Susan_Johnston_Grant_Figures/Males_figure.pdf", height = 15, width = 15)
system("open data//Susan_Johnston_Grant_Figures/Males_figure.pdf")

cell_numbers <- lapply(plots, function(x)(return(x[[3]]))) %>% 
  bind_rows()
cell_numbers %>% 
  filter(Freq != 0) %>% 
  write.csv("data/Susan_Johnston_Grant_Figures/cell_numbers.csv")
