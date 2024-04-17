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

germ_somatic_files <- list.files("data/mut_rate/germ_somatic/", 
                                 full.names = T, pattern = "germ")
sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")


x <- read.table(germ_somatic_files[1])


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
    mutate(sample = gsub("_", "", sample) %>% toupper()) %>% 
    merge(si) %>% 
    mutate(newname = paste(sex, ed, sep = " "), 
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
depths = c(5,6,7,8,9,10,12,15,20,25,30)
depths = 30
check_depths = do.call(rbind, lapply(germ_somatic_files, nest_dsf, dps = depths))

mut_level_plot_one <- function(dp, cd){
  cd <- filter(cd, somatic_depth == dp)
  comps <- check_depths$newname %>% unique() %>% combn(2) %>% t() %>% 
    split(.,seq(nrow(.)))
  
  
  plot <- cd %>% 
    filter(mut_level > 0) %>% 
    ggplot(aes(x = newname, y = logml)) + geom_boxplot() + 
    stat_compare_means(comparisons = comps, method = 'wilcox.test') + 
    geom_jitter(size = 0.01) + 
    ggtitle(paste("depth = ", dp))
  return(plot)
}

plots <- lapply(depths, mut_level_plot, cd = check_depths)
names(plots) <- depths
arranged_plots <- ggarrange(plotlist = plots)
ggsave("plots/arranged.pdf", arranged_plots, height = 40, width = 40)

mut_level_plot_bar <- function(dp, cd){
  cd <- filter(cd, somatic_depth == dp)
  plot <- cd %>% 
    group_by(sex, ed, sample, mutant, newname) %>% 
    summarise(numbers = n()) %>% 
    ggplot(aes(x = newname, y = numbers, fill = sample, colour = mutant)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    ggtitle(paste("depth = ", dp))
  return(plot)
}

plots_bar <- lapply(depths, mut_level_plot_bar, cd = check_depths)
names(plots_bar) <- depths
arranged_plots <- ggarrange(plotlist = plots_bar)
ggsave("plots/arranged_bar.pdf", arranged_plots, height = 40, width = 40)





