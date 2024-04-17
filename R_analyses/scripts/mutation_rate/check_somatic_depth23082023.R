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

metadata_get <- function(x){
  load(x)
  keep_cols <- 
    c("nUMI", "nGene", "cells", "log10GenesPerUMI", "mitoRatio", 
      "nCount_SCT", "nFeature_SCT", 
      "Phase", "S.Score", "G2M.Score", "seurat_clusters", "integrated_snn_res.3")
  return(seurat_marker@meta.data[,keep_cols])
}

seurat_metadata <- do.call(rbind, lapply(RDatas, metadata_get)) 

germ_somatic_files <- list.files("data/mut_rate/germ_somatic2230823//", 
                                 full.names = T, pattern = "germ")
sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"


x <- read.table(germ_somatic_files[1])
dp=10
si=sample_info
smd=seurat_metadata

depth_summary_func <- function(dp, x, si, smd=seurat_metadata){
  gs <- x %>% 
    filter(genotype.other == 'het' | genotype.somatic != 'het') %>%
    filter(depth.somatic >=dp) %>% 
    group_by(Barcode, sample, sctype_labels) %>% 
    reframe(nmuts = length(which(mutated == 'mutation')),
            nnone = length(which(mutated == "none")),
            mut_level = nmuts/(nnone+nmuts), 
            somatic_depth = dp, 
            logml = log(mut_level + 0.001)) %>% 
    merge(si) %>% 
    mutate(comb_name = paste(sex, stage, sep = " "), 
           mutant = mut_level > 0) %>% 
    merge(.,smd, by.x = "Barcode", by.y = "cells") %>% 
    filter((sex == "F" & ed == "ED24" & 'integrated_snn_res.3' %in% c(26,22,23)) == F)
  return(gs)
}


nest_dsf <- function(file, dps = c(15), si = saple_info){
  gs_tmp <- read.table(file)
  tmp <- do.call(rbind, lapply(dps, depth_summary_func, x = gs_tmp, 
                               si = sample_info))
  return(tmp)
}
depths = c(5,7,10)
check_depths = do.call(rbind, lapply(germ_somatic_files, nest_dsf, dps = depths))

mut_level_plot <- function(dp, cd){
  cd <- filter(cd, somatic_depth == dp)
  comps <- check_depths$comb_name %>% unique() %>% combn(2) %>% t() %>% 
    split(.,seq(nrow(.)))
  
  
  plot <- cd %>% 
    filter(mut_level > 0) %>% 
    ggplot(aes(x = comb_name, y = logml)) + geom_boxplot() + 
    stat_compare_means(comparisons = comps, method = 'wilcox.test') + 
    geom_jitter(size = 0.01) + 
    ggtitle(paste("depth = ", dp)) + 
    facet_grid(. ~ species)
  return(plot)
}

plots <- lapply(depths, mut_level_plot, cd = check_depths)
names(plots) <- depths
arranged_plots <- ggarrange(plotlist = plots)
ggsave("plots/arranged.pdf", arranged_plots, height = 40, width = 40)

mut_level_plot_bar <- function(dp, cd){
  cd <- filter(cd, somatic_depth == dp)
  plot <- cd %>% 
    group_by(sex, ed, sample, mutant, comb_name, species) %>% 
    summarise(numbers = n()) %>% 
    ggplot(aes(x = comb_name, y = numbers, fill = sample, colour = mutant)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    ggtitle(paste("depth = ", dp)) +
    facet_grid(. ~ species)
  return(plot)
}

plots_bar <- lapply(depths, mut_level_plot_bar, cd = check_depths)
names(plots_bar) <- depths
arranged_plots <- ggarrange(plotlist = plots_bar)
ggsave("plots/arranged_bar.pdf", arranged_plots, height = 40, width = 40)






mut_level_plot_bar_prop <- function(dp, cd){
  cd <- filter(cd, somatic_depth == dp)
  plot <- cd %>% 
    group_by(sex, ed, sample, comb_name, species) %>% 
    summarise(muts = length(which(mutant == TRUE)), 
              none = length(which(mutant == FALSE))) %>% 
    mutate(prop_mut_cells = muts/(muts + none)) %>% 
    ggplot(aes(x = comb_name, y = prop_mut_cells, fill = sample)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    ggtitle(paste("depth = ", dp)) +
    facet_grid(. ~ species)
  return(plot)
}

plots_bar_prop <- lapply(depths, mut_level_plot_bar_prop, cd = check_depths)
names(plots_bar_prop) <- depths
arranged_plots <- ggarrange(plotlist = plots_bar_prop)
ggsave("plots/arranged_bar_prop.pdf", arranged_plots, height = 40, width = 40)

plots <- list()
plots[[1]] <- cd %>% 
  group_by(sex, ed, sample, comb_name, species) %>% 
  summarise(muts = length(which(mutant == TRUE)), 
            none = length(which(mutant == FALSE))) %>% 
  mutate(prop_mut_cells = muts/(muts + none)) %>% 
  ggplot(aes(x = comb_name, y = prop_mut_cells, fill = sample)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  ggtitle(paste("depth = ", dp)) +
  facet_grid(. ~ species) + theme(legend.position = "none")
plots[[2]] <-cd %>% 
  group_by(sex, ed, sample, comb_name, species) %>% 
  summarise(muts = length(which(mutant == TRUE)), 
            none = length(which(mutant == FALSE))) %>% 
  mutate(prop_mut_cells = muts/(muts + none)) %>% 
  ggplot(aes(x = comb_name, y = prop_mut_cells, fill = comb_name)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps[2:5], method = 'wilcox.test') + 
  ggtitle(paste("depth = ", dp)) +
  facet_grid(. ~ species) + theme(legend.position = "none")
plots[[3]] <-cd %>% 
  group_by(sex, stage, sample, comb_name, species) %>% 
  #summarise(muts = length(which(mutant == TRUE)), 
  #          none = length(which(mutant == FALSE))) %>% 
  #mutate(mut_level = muts/(muts + none)) %>% 
  filter(mut_level > 0) %>% 
  ggplot(aes(x = comb_name, y = logml, fill = sex)) + 
  geom_boxplot()  +
  stat_compare_means(comparisons = comps[2:5], method = 'wilcox.test') + 

  #geom_jitter(size = 0) +
  ggtitle(paste("depth = ", dp)) +
  facet_grid(. ~ species)

arranged <- ggarrange(plotlist = plots, nrow = 3)
ggsave("plots/mutation_levels.pdf", arranged)
RDatas <- list.files("data/seurat_RData/", full.names = T)

plot_PCA_and_markers <- function(x){
  load(x)
  
  si_tmp <- filter(sample_info, sample %in% unique(seurat_marker$sample))[1,]
  save_name <- paste(si_tmp$species, si_tmp$sex, si_tmp$stage, sep = "_")
  plots <- list()
  
  DefaultAssay(seurat_marker) <- "RNA"
  plots[[1]] <- FeaturePlot(seurat_marker, features = "DAZL") + ggtitle("DAZL germ")
  plots[[2]] <- FeaturePlot(seurat_marker, features = "DDX4")+ ggtitle("DDX4 germ")
  plots[[3]] <- FeaturePlot(seurat_marker, features = "AMH") + ggtitle("AMH sertoli")
  plots[[4]] <- FeaturePlot(seurat_marker, features = "SOX18") + ggtitle("SOX18 epithelial")
  plots[[5]] <- FeaturePlot(seurat_marker, features = "FOXL2") + ggtitle("FOXL2 granulosa")
  plots[[6]] <- FeaturePlot(seurat_marker, features = "TAGLN3") + ggtitle("TAGLN3 neuronal")
  plots[[7]] <- FeaturePlot(seurat_marker, features = "POSTN") + ggtitle("POSTN interstitial")
  
  plots[[8]] <- DimPlot(seurat_marker, group.by = 'sctype_labels') + ggtitle(save_name)
  #plots[[9]] <- DimPlot(seurat_marker, group.by = 'scina_labels') + ggtitle(save_name)
  
  arranged <- ggarrange(plotlist = plots)

  ggsave(paste("plots/UMAP_", save_name, ".pdf", sep = ""), arranged, height = 15, width = 15)
}

lapply(RDatas, plot_PCA_and_markers)
