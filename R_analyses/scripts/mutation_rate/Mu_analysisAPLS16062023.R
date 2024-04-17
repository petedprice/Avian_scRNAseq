library("optparse")
library(SCINA)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)


#SC_TYPE FUNCTIONS
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('Sample Name', "species", 'ref', "Sex", "ed")
sample_info$files <- gsub("APLS", "apls_", sample_info$`Sample Name`)
sample_info <- filter(sample_info, species == 'duck')
load("data/mut_rate/mut_seurat.R")

update_rownames <- function(x){
  rownames(x@meta.data) <- x@meta.data$cells
  return(x)
}

mut_seurat <- lapply(mut_seurat, update_rownames)

mut_metadata <- function(x){
  temp <- x@meta.data
  temp <- filter(temp, sctype_labels == 'germ')
  return(temp)
}
 
mut_seurat_germ <- do.call(rbind, lapply(mut_seurat, mut_metadata))
mut_seurat_germ <- merge(mut_seurat_germ, sample_info[,c(4,5,6)], by.x = 'sample', by.y = 'files') 
mut_seurat_germ$newname <- paste(mut_seurat_germ$Sex, mut_seurat_germ$ed, sep = " ")

mut_seurat_germ %>% ggplot(aes(x = ed, color = Sex, y = mut_level)) + geom_boxplot() + 
  stat_compare_means(method = 't.test')
comps <- list()
comps[[1]] <- c("F ED17", "F ED24")
comps[[2]] <- c("M ED17", "M ED24")
comps[[3]] <- c("F ED17", "M ED17")
comps[[4]] <- c("F ED24", "M ED24")

mut_seurat_germ %>% ggplot(aes(x = newname, y = mut_level, fill = sctype_labels)) + geom_boxplot() + 
  stat_compare_means(comparisons = comps)

lm(mut_level ~ Sex + ed, mut_seurat_germ) %>% 
  summary()

mut_levelz <- mut_seurat_germ %>% 
  group_by(sctype_labels, newname) %>% 
  summarise(mutlevel = mean(mut_level, na.rm = T),
            mutvar = var(mut_level, na.rm = T), 
            se = sd(mut_level, na.rm = T)/sqrt(length(n)),
            sd = sd(mut_level, na.rm = T),
            n = n())

mut_levelz %>% 
  ggplot(aes(x = newname, y = mutlevel, fill = sctype_labels)) + geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin=mutlevel, ymax=mutlevel+mutvar), width=.2,
                position=position_dodge(.9))
mut_seurat_germ %>% filter(ed == 'ED24') %>% 
  ggplot(aes(x = nnone, y = nGene, color = Sex)) + geom_point()

mut_seurat_germ %>% filter(ed == 'ED24') %>% 
  ggplot(aes(x = nUMI, color = Sex)) + geom_density()

mut_seurat_germ %>% group_by(ed, Sex) %>% 
  summarise(n = n())


d <- lapply(mut_seurat, function(x)(ggplot(x@meta.data, aes(x = nGene)) + geom_density() + ggtitle(x@meta.data$sample)))
ggarrange(plotlist = d)


ggplot(filter(mut_seurat$apls_34@meta.data), aes(x = nGene, color = sctype_labels)) + geom_density()

Fe(mut_seurat$apls_34, group.by = 'sctype_labels')
DimPlot(mut_seurat$apls_3, label = T)
DimPlot(mut_seurat$apls_3, group.by = 'sctype_labels')
FeaturePlot(mut_seurat$apls_3, features = c("nUMI"))

