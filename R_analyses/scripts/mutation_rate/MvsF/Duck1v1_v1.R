# MALE VS FEMALE 
###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)

obj_list <- list()
files = c("apls_29", "apls_30")
for (file in files){
  print(file)
  
  sample <- file
  
  seurat_data <- Read10X(data.dir = paste("data/cellranger/", file, "/outs/filtered_feature_bc_matrix/", sep = ""))
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200, ### FILTERING CELLS WITH LESS THAN 200 GENES
                                   project = sample, 
                                   min.cells = 3) #Might be able to include min.cells = 3 for keeping a gene rather than doing this later
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample)
  seurat_obj@meta.data$sample <- sample
  seurat_obj@meta.data$treatment <- substr(sample, 1, 2)
  obj_list <- c(obj_list, seurat_obj)
}
merged_seurat <- reduce(obj_list, merge)
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA) # Calculating the number of features per UMI 


#CALCULATING PROPORTION OF READS MAPPING TO MITOCHONDRIAL GENOME (x100 as originally a percentage)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "mt")/ 100 
metadata <- merged_seurat@meta.data # Seperately creating metadata dataframe to save metrics etc without risking affecting seurat object
metadata$cells <- rownames(metadata) #add cell IDs to metadat
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
merged_seurat@meta.data <- metadata #Save the more complete metadat to the seurat object

#FILTERING DATA
filtered_seurat <- subset(x = merged_seurat, 
                          (log10GenesPerUMI > 0.8) & # Can be dying cells or simple cell types such as blood cells
                            (mitoRatio < 0.05))
metadata_clean <- filtered_seurat@meta.data

filtered_seurat <- SCTransform(filtered_seurat, vars.to.regress = 
                                 c("mitoRatio","nUMI"))

split_seurat <- SplitObject(filtered_seurat, split.by = 'sample')

seurat_SCT_normalised <- sapply(split_seurat, RunPCA)
seurat_SCT_normalised <- sapply(seurat_SCT_normalised, RunTSNE, reduction  = 'pca')

seurat_SCT_normalised <- sapply(seurat_SCT_normalised, RunUMAP, 
                                dims = 1:30,
                                reduction = "pca", 
                                resolution = 0.2)
seurat_SCT_normalised <- sapply(seurat_SCT_normalised, FindNeighbors,  
                                dims = 1:30, verbose = FALSE)

seurat_SCT_normalised <- sapply(seurat_SCT_normalised, FindClusters, 
                                verbose = FALSE, resolution = 0.4)

sapply(seurat_SCT_normalised, DimPlot)

DimPlot(seurat_SCT_normalised$apls_29, group.by = 'ident', label = T)
DimPlot(seurat_SCT_normalised$apls_30, group.by = 'ident', label = T)

FeaturePlot(seurat_SCT_normalised$apls_29, features = 'DAZL')
FeaturePlot(seurat_SCT_normalised$apls_30, features = 'DAZL')

germ_markers <- c("DAZL", "DDX4", "CLDN1", "DND1", "LOC101795050", "CARHSP1", "Vasa")
FeaturePlot(seurat_SCT_normalised$apls_29, germ_markers)
FeaturePlot(seurat_SCT_normalised$apls_30, germ_markers)

somatic_markers <- c("POSTN", "LOC101798048", "SOX9", "DMRT1")
FeaturePlot(seurat_SCT_normalised$apls_29, somatic_markers)
FeaturePlot(seurat_SCT_normalised$apls_30, somatic_markers)

## SAMPLE 29 germ cells are 2
## SAMPLE 30 germ cells are 3,5,7

## SAMPLE 29 SOMATIC cells are 1,3,7
## SAMPLE 30 SOMATIC cells are 6,2

seurat_SCT_normalised$apls_29$type <- seurat_SCT_normalised$apls_29$SCT_snn_res.0.4 %>% as.numeric()
#seurat_SCT_normalised$apls_29$type[seurat_SCT_normalised$apls_29$SCT_snn_res.0.4 %in% c(2)] <- "germ"
seurat_SCT_normalised$apls_29$type[seurat_SCT_normalised$apls_29$SCT_snn_res.0.4 %in% c(1,3,7)] <- "somatic"
#seurat_SCT_normalised$apls_29$type[seurat_SCT_normalised$apls_29$SCT_snn_res.0.4 %in% c(1,2,3,7) == F] <- "other"


seurat_SCT_normalised$apls_30$type <- seurat_SCT_normalised$apls_30$SCT_snn_res.0.4 %>% as.numeric()
#seurat_SCT_normalised$apls_30$type[seurat_SCT_normalised$apls_30$SCT_snn_res.0.4 %in% c(3,5,7)] <- "germ"
seurat_SCT_normalised$apls_30$type[seurat_SCT_normalised$apls_30$SCT_snn_res.0.4 %in% c(6,2)] <- "somatic"
#seurat_SCT_normalised$apls_30$type[seurat_SCT_normalised$apls_30$SCT_snn_res.0.4 %in% c(2,3,5,6,7) == F] <- "other"


mut_rate29 <- read.table("data/mut_rate/apls_29_mutation_data2.txt.gz", sep = " ", header = F) %>% 
  select(-c(V1))
mut_rate30 <- read.table("data/mut_rate/apls_30_mutation_data2.txt.gz", sep = " ", header = F) %>% 
  select(-c(V1))

colnames(mut_rate29) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")
colnames(mut_rate30) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")
mut_rate29$Barcode <- paste("apls_29_", mut_rate29$Barcode, sep = "")
mut_rate30$Barcode <- paste("apls_30_", mut_rate30$Barcode, sep = "")

mut_rate29_sc <- mut_rate29 %>% merge(seurat_SCT_normalised$apls_29@meta.data[
  is.na(seurat_SCT_normalised$apls_29$type) == F,c('cells', 'type')], 
                                         by.x = 'Barcode', by.y = "cells")
#som_number <- 2500000
#germ_number <- 6000000

mut_rate30_sc <- mut_rate30 %>% merge(seurat_SCT_normalised$apls_30@meta.data[
  is.na(seurat_SCT_normalised$apls_30$type) == F,c('cells', 'type')], 
                                   by.x = "Barcode", by.y = "cells")

#s29 <- sample(rownames(mut_rate29_sc)[mut_rate29_sc$type == 'somatic'], som_number)
#s30 <- sample(rownames(mut_rate30_sc)[mut_rate30_sc$type == 'somatic'], som_number)
#g29 <- sample(rownames(mut_rate29_sc)[mut_rate29_sc$type == 'germ'], germ_number)
#g30 <- sample(rownames(mut_rate30_sc)[mut_rate30_sc$type == 'germ'], germ_number)


mut_rate29_sc$sample <- "apls29"
mut_rate30_sc$sample <- "apls30"

merged_sc_mut_save <- rbind(mut_rate30_sc, mut_rate29_sc)
#merged_sc_mut <- merged_sc_mut_save %>% group_by(sample, type) %>% 
#  slice_sample(n = 2500000)
merged_sc_mut <- merged_sc_mut_save 

merged_sc_mut %>% ggplot(aes(x = type, fill = sample)) + geom_bar(position = 'dodge')

somatic_snps <- filter(merged_sc_mut, type == "somatic") %>% 
  group_by(pos, CHROM, sample) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  filter(altn > 1 | refn > 1) %>% 
  #filter(altn == depth | refn == depth) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref", 
                              refn > 0 & altn > 0 ~ "het")) %>% 
  filter(genotype != 'het')

#??????????????????????????????????
germ_snps <- filter(merged_sc_mut, type != "somatic") %>% 
  group_by(pos, Barcode, sample, type) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  #filter(altn > 1 & refn > 1) %>% 
  filter(altn == depth | refn == depth) %>% 
  filter(pos %in% somatic_snps$pos) %>% 
  #mutate(genotype = 'het') %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref")) %>% 
  merge(somatic_snps, by = c('pos', 'sample'), suffixes = c( '.germ','.somatic')) %>% 
  mutate(mutated = case_when(genotype.germ != genotype.somatic ~ "mutation", 
                             genotype.germ == genotype.somatic ~ "none")) %>% 
  filter(is.na(mutated) == F) %>% 
  group_by(pos) %>%
  filter(n() == 1)


gs <- germ_snps %>% 
  filter(depth.germ > 0) %>% 
  group_by(Barcode, sample, type) %>% 
  reframe(nmuts = sum(which(mutated == 'mutation')),
          nnone = sum(which(mutated == "none")))
gs$mut_level <- gs$nmuts/(gs$nnone+gs$nmuts)


gs %>% ggboxplot(x = 'type', y = 'mut_level', color = 'type', add = 'jitter') + 
  stat_compare_means(method = 't.test')
plots <- list()
plots[[1]] <- gs %>% ggplot(aes(y = mut_level, x = type, fill = type))  + 
  geom_boxplot() 
plots[[2]] <- DimPlot(seurat_SCT_normalised$apls_29, group.by = 'type', label = T)  

ggarrange(plotlist = plots)

gs %>% ggplot(aes(x = mut_level, colour = type)) + facet_grid(~sample) + 
  geom_density()

mut_levelz <- gs %>% 
  group_by(type, sample) %>% 
  summarise(mutlevel = mean(mut_level),
            mutvar = var(mut_level), 
            se = sd(mut_level)/sqrt(length(n)),
            sd = sd(mut_level),
            n = n())

mut_levelz %>% 
  ggplot(aes(x = type, y = mutlevel, fill = type)) + geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin=mutlevel, ymax=mutlevel+sd), width=.2,
                position=position_dodge(.9)) 

