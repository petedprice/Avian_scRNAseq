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

mut_files <- list.files("data/mut_rate/mut_metadata28092023/mutant_metadata28092023/", full.names = T)
sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")


read_mut_csv <-function(x){
  tmp <- read.csv(x)
  keep_cols <- c("nUMI", "nGene", "cells", "log10GenesPerUMI", "mitoRatio", 
                 "nCount_SCT", "nFeature_SCT", "sample", "nmuts", "nnone", "mut_level", 
                 "Phase", "S.Score", "G2M.Score", "seurat_clusters")
  return(tmp[,keep_cols])
}

mut_data <- do.call(rbind, lapply(mut_files, read_mut_csv)) %>% 
  merge(sample_info, by = 'sample') %>% 
  mutate(logml = log(mut_level + 0.00001), 
         nm = nnone + nmuts, 
         newname = paste(sex, ed, sep = " "), 
         mutant = mut_level > 0)
  #filter(sex != "F" | ed != "ED24" | seurat_clusters != 3)
  #filter(logml < mean(logml) + (3 * sd(logml)) & logml > mean(logml) - (3 * sd(logml)))



update_obj <- function(obj, mut_data){
  tmp_md <- filter(mut_data, cells %in% colnames(obj))
  obj@meta.data <- 
    merge(obj@meta.data,
          mut_data[,c('cells', "nmuts", "nnone", "mut_level", "logml", "mutant")], all.x = T)
  obj@meta.data <- obj@meta.data %>% 
    mutate(logumi = log(nUMI), 
           dying = (log10GenesPerUMI < 0.8), 
           mutated = mut_level > 0, 
           newname = tmp_md$newname[1])
  rownames(obj@meta.data) <- obj@meta.data$cells
  
  
  return(obj)
}

obj_list <- list()
load("data/seurat_RData/duck_F_ED17_marker_seurat.RData")
obj_list[['seurat_F17']] <- seurat_marker %>% update_obj(., mut_data)

load("data/seurat_RData/duck_F_ED24_marker_seurat.RData")
obj_list[['seurat_F24']] <- seurat_marker %>% update_obj(., mut_data)

load("data/seurat_RData/duck_M_ED17_marker_seurat.RData")
obj_list[['seurat_M17']] <- seurat_marker %>% update_obj(., mut_data)

load("data/seurat_RData/duck_M_ED24_marker_seurat.RData")
obj_list[['seurat_M24']] <- seurat_marker %>% update_obj(., mut_data)


comps <- list()
comps[[1]] <- c("F ED17", "F ED24")
comps[[2]] <- c("M ED17", "M ED24")
comps[[3]] <- c("F ED17", "M ED17")
comps[[4]] <- c("F ED24", "M ED24")


mut_data %>% 
  filter(mut_level > 0) %>% 
  ggplot(aes(x = newname, y = logml)) + geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test') + 
  geom_jitter(size = 0.01)

mut_data %>% 
  filter(mut_level > 0) %>% 
  lmerTest::lmer(logml ~ sex * ed + (1|nGene) + 
                  (1|G2M.Score) + 
                 (1|S.Score), .) %>% summary()

mut_data %>% 
  filter(mut_level > 0) %>% 
  ggplot(aes(x = newname, y = logml, fill = sample)) + geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test') + 
  geom_jitter(size = 0.01)

mut_data %>% 
  filter(mut_level > 0) %>% 
  lmerTest::lmer(logml ~ sex * ed + (1|nGene) + (1|sample) +
                    (1|G2M.Score) + 
                   (1|S.Score), .) %>% summary()


mut_data %>% 
  group_by(sex, ed, sample, mutant, newname) %>% 
  summarise(numbers = n()) %>% 
  ggplot(aes(x = newname, y = numbers, fill = sample, colour = mutant)) + 
  geom_bar(stat = 'identity', position = 'dodge')

mut_data %>% 
  glmer(mutant ~ sex * ed + (1|nGene) + (1|sample) +
          (1|mitoRatio) + (1|G2M.Score) + 
          (1|S.Score), ., family = 'binomial') %>% 
  summary()



mut_data_summary <- mut_data %>% 
  dplyr::count(mutant, ed, sex, sample) %>% 
  spread(mutant, n)

glmer(cbind(mut_data_summary$MUTANT, 
            mut_data_summary$MUTANT + mut_data_summary$NOPE) ~ 
        sex * ed + (1|sample), data = mut_data_summary, family = "binomial") %>% 
  summary()

UMAPS <- lapply(obj_list$seurat_F24, function(x)(DimPlot(x, group.by = 'sctype_labels', label = T) + ggtitle(x$newname[1])))
ggarrange(plotlist = UMAPS)

DF_assay <- function(x){
  DefaultAssay(x) <- "RNA"
  return(x)
}
obj_list <- lapply(obj_list$seurat_F24, DF_assay)
feature_list <- lapply(obj_list$seurat_F24, function(x)(
    FeaturePlot(x,features = c("DAZL", "nGene", "logml", "nUMI")) + ggtitle(x$newname[1])))

ggarrange(plotlist = feature_list)
