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

#AverageExpression() calcaultes average expression of a classification of celltype
exp_avs <- lapply(seur_objs, AverageExpression, group.by = 'sctype_labels', assays = "RNA")

get_gene_names <- function(x){
  return(rownames(x[[1]]))
}

gene_names <- lapply(exp_avs, get_gene_names)
genes <- unlist(gene_names) %>% table() %>% as.data.frame() %>% filter(Freq > 11)
colnames(genes) <- c("Gene", "Freq")
gene <- genes$Gene[1]

samp_names <- names(exp_avs)

add_sample_info <- function(name, exp_avs){
  tmp <- exp_avs[[name]][[1]] %>% as.data.frame()
  tmp <- tmp %>% mutate(gene = rownames(tmp))
  tmp <- tmp %>% pivot_longer(colnames(tmp)[-ncol(tmp)], names_to = "cell_type", 
                              values_to = "expression")
  details <- str_split(name, "_", simplify = T)
  colnames(details) <- c("species", "sex", "stage")
  tmp <- cbind(tmp, details)
  return(tmp)
}

compiled_data <- do.call(rbind, lapply(samp_names, add_sample_info, exp_avs))

gene_variances <- function(goi, x){
  data <- subset(x, gene == goi & cell_type != "Unkown")
  variance <- var(data$expression)
  return(variance)
}

variances <- sapply(genes$Gene, gene_variances, x = compiled_data)
names(variances) <- genes$Gene

gene_means <- function(goi, x){
  data <- subset(x, gene == goi & cell_type != "Unkown")
  variance <- mean(data$expression)
  return(variance)
}

means <- sapply(genes$Gene, gene_means, x = compiled_data)
names(means) <- genes$Gene

varmeans <- data.frame(mu = means, sigma = variances)
varmeans$scaled_var <- varmeans$sigma/varmeans$mu


tmp <- filter(varmeans, mu > 1)
a <- tmp[order(tmp$scaled_var, decreasing = F),] 
View(a)

b <- compiled_data[compiled_data$gene == 'HECA',]
View(b)

variances[order(variances)] %>% names()







