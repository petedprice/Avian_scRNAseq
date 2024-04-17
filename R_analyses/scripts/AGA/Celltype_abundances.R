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

sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"

load_func <- function(x){
  load(x)
  return(seurat_marker)
}

RDatas <- list.files("data/seurat_RData/", full.names = T)
seur_objs <- lapply(RDatas, load_func)
length(seur_objs)
names(seur_objs) <- str_split(RDatas, "/", simplify = T)[,4] %>% 
  gsub("_marker_seurat.RData", "", .)


cell_type_abund_func <- function(seur_obj){
  abundances <- table(seur_obj@meta.data$sctype_labels)
  return(abundances)
  
}

abundances_all_species <- do.call(bind_rows, lapply(seur_objs, cell_type_abund_func)) %>% 
  as.data.frame()
abundances_all_species <- abundances_all_species[,-6]
rownames(abundances_all_species) <- names(seur_objs)
sums <- rowSums(abundances_all_species, na.rm = T)

aas_prop <- apply(abundances_all_species, 2, function(x)(return(x/sums)))
aas_prop

      