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
library(future)
library(future.apply)
plan("multicore", workers = 2)
#options(future.globals.maxSize = 8000 * 1024^5)

#SC_TYPE FUNCTIONS
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


#sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
sample_info <- read.csv("/home/bop20pp/software/Avian_scRNAseq/CL_analyses/nextflow/bioworkbook/metadata_APLS.csv", header = F) 

colnames(sample_info) <- c('Sample Name', "species", 'ref', "Sex", "ed")
sample_info$files <- gsub("APLS", "apls_", sample_info$`Sample Name`)
#sample_info <- filter(sample_info, species == 'duck' & 
#                        files %in% c('apls_1', 'apls_2', 'apls_12', 'apls_16')
#)
obj_list <- list()
for (file in sample_info$files){
  print(file)
  
  sample <- file
  
  #seurat_data <- Read10X(data.dir = paste("data/cellranger/", file, "/outs/filtered_feature_bc_matrix/", sep = ""))
  seurat_data <- Read10X(data.dir = paste("/fastdata/bop20pp/Avian_scRNAseq/cellranger/", file, "/outs/filtered_feature_bc_matrix/", sep = ""))
  
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
split_seurat <- SplitObject(filtered_seurat, split.by = 'sample')
seurat_SCT_normalised <- future_lapply(split_seurat, SCTransform, vars.to.regress = 
                                  c("mitoRatio","nUMI"))

seurat_SCT_normalised <- future_lapply(seurat_SCT_normalised, RunPCA)
seurat_SCT_normalised <- future_lapply(seurat_SCT_normalised, RunTSNE, reduction  = 'pca')

seurat_SCT_normalised <- future_lapply(seurat_SCT_normalised, RunUMAP, 
                                dims = 1:30,
                                reduction = "pca", 
                                resolution = 0.2)
seurat_SCT_normalised <- future_lapply(seurat_SCT_normalised, FindNeighbors,  
                                dims = 1:30, verbose = FALSE)

seurat_SCT_normalised <- future_lapply(seurat_SCT_normalised, FindClusters, 
                                verbose = FALSE, resolution = 0.4)
markerlist <- list()
markerlist[['interstitial']] <- c("POSTN", "DCN")
markerlist[['germ']] <- c("DAZL","DDX4")
markerlist[['sertoli']] <- c("AMH")
markerlist[['epethlial']] <- c("SOX18", "CDH5")
markerlist[['neuronal']] <- c("TAGLN3")
markerlist[['GRANULOSA']] <- c("FOXL2", "FGFR3")
gs_list <- list()
gs_list$gs_positive <- markerlist


marker_function <- function(seur_obj, gs_list){
  es.max = sctype_score(scRNAseqData = seur_obj[["SCT"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(seur_obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seur_obj@meta.data[seur_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seur_obj@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  seur_obj@meta.data$sctype_labels = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seur_obj@meta.data$sctype_labels[seur_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  return(seur_obj)
}

seurat_SCT_normalised <- sapply(seurat_SCT_normalised, marker_function, gs_list = gs_list)


plotlist = lapply(c(1:nrow(sample_info)), function(x)(return(
  DimPlot(seurat_SCT_normalised[[sample_info$files[x]]], group.by = 'sctype_labels', label = T) + 
    ggtitle(paste(sample_info$files[x], sample_info$ed[x], sample_info$Sex[x], sep = "_"))
)))

for (s in unique(sample_info$Sex)){
  for (day in unique(sample_info$ed)){
    ss <- filter(sample_info, Sex == s & ed == day)
    plotlist = lapply(c(1:nrow(ss)), function(x)(return(
      DimPlot(seurat_SCT_normalised[[ss$files[x]]], group.by = 'sctype_labels', label = T) + 
        ggtitle(paste(ss$files[x], ss$ed[x], ss$Sex[x], sep = "_"))
    )))
    plots <- ggarrange(plotlist = plotlist)
    ggsave(paste("APLS_UMAP", day, s, ".pdf", sep = ""), plots, height = 12, width = 12)
  }
}


plots <- ggarrange(plotlist = plotlist)
ggsave("APLS_UMAP.pdf", plots, height = 40, width = 40)

mutation_read <- function(sample){
  #mut_rate <- read.table(paste("data/mut_rate/", sample, "_mutation_data.txt.gz", sep = ""), sep = " ", header = F) %>% 
  mut_rate <- read.table(paste("/fastdata/bop20pp/Avian_scRNAseq/nextflow/mut_compiled/send/", sample, "_mutation_data.txt.gz", sep = ""), sep = " ", header = F) 
    #select(-c(V1))
  colnames(mut_rate) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")
  mut_rate$Barcode <- paste(sample, mut_rate$Barcode, sep = "_")
  
  mut_rate_sc <- mut_rate %>% merge(seurat_SCT_normalised[[sample]]@meta.data[
    is.na(seurat_SCT_normalised[[sample]]$sctype_labels) == F,c('cells', 'sctype_labels')], 
    by.x = 'Barcode', by.y = "cells")
  mut_rate_sc$sample <- sample
  #write.csv(mut_rate_sc,  file=gzfile(paste(sample, "merged_sc_mut.csv.gz", sep = "_")))
  
  return(mut_rate_sc)
}
c=0
for (f in sample_info$files){
  c=c+1
  print(f)
  if (c==1){
    merged_temp <- mutation_read(f)
  } else {
    temp <- mutation_read(f)
    merged_temp <- rbind(merged_temp, temp)
  }
}


#merged_sc_mut <- do.call(rbind, merged_temp)
merged_sc_mut <- merged_temp
#write.csv(merged_sc_mut,  file=gzfile("merged_sc_mut.csv.gz"))

somatic_snps <- filter(merged_sc_mut, sctype_labels != "germ") %>% 
  group_by(pos, CHROM, sample) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref", 
                              refn > 0 & altn > 0 ~ "het")) 

write.csv(somatic_snps,  file=gzfile("/fastdata/bop20pp/Avian_scRNAseq/R_analyses/mut_id/somatic_snps.csv.gz"))

somatic_snps_filt <- somatic_snps %>% 
  filter(altn > 10 | refn > 10) %>% 
  filter(altn == depth | refn == depth) %>% 
  filter(genotype != 'het')
  
#??????????????????????????????????
germ_snps <- filter(merged_sc_mut, sctype_labels == "germ") %>% 
  group_by(pos, CHROM, Barcode, sample, sctype_labels) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  filter(pos %in% somatic_snps$pos) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref",
                              (refn > 0 & altn > 0) ~"het"))

write.csv(germ_snps,  file=gzfile(
  "/fastdata/bop20pp/Avian_scRNAseq/R_analyses/mut_id/germ_snps.csv.gz"))

germ_somatic <-  germ_snps %>% 
  merge(somatic_snps, by = c('pos', 'sample', 'CHROM'), suffixes = c( '.other','.somatic')) %>% 
  mutate(mutated = case_when(genotype.other != genotype.somatic ~ "mutation", 
                             genotype.other == genotype.somatic ~ "none"))

germ_somatic$mutated[germ_somatic$genotype.other == 'het' & 
                       germ_somatic$genotype.somatic != 'het'] <- "mutation"


gs <- germ_somatic %>% 
 filter(genotype.other == 'het' | genotype.somatic != 'het') %>% 
  group_by(Barcode, sample, sctype_labels) %>% 
  reframe(nmuts = sum(which(mutated == 'mutation')),
          nnone = sum(which(mutated == "none")))
gs$mut_level <- gs$nmuts/(gs$nnone+gs$nmuts)


write.table(gs, "//fastdata/bop20pp/Avian_scRNAseq/R_analyses/mut_id/gs.txt")



gs <-read.table("data/mut_rate/gs.txt") 
gs <- merge(gs, sample_info[,c(4,5,6)], by.x = 'sample', by.y = 'files') %>% 
  filter(log(nnone + nmuts) > 10)
gs$newname <- paste(gs$Sex, gs$ed, sep = " ")

comps <- list()
comps[[1]] <- c("F ED17", "F ED24")
comps[[2]] <- c("M ED17", "M ED24")
comps[[3]] <- c("F ED17", "M ED17")
comps[[4]] <- c("F ED24", "M ED24")

gs %>% ggplot(aes(x = newname, y = mut_level, fill = sctype_labels)) + geom_boxplot() + 
  stat_compare_means(comparisons = comps)
gs %>% ggplot(aes(x = mut_level, color = sample)) + geom_density()

mut_levelz <- gs %>% 
  group_by(sctype_labels, newname) %>% 
  summarise(mutlevel = mean(mut_level),
            mutvar = var(mut_level), 
            se = sd(mut_level)/sqrt(length(n)),
            sd = sd(mut_level),
            n = n())

mut_levelz %>% 
  ggplot(aes(x = newname, y = mutlevel, fill = sctype_labels)) + geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin=mutlevel, ymax=mutlevel+sd), width=.2,
                position=position_dodge(.9)) +
  stat_compare_means(method = 't.test')

gs %>% 
  filter(sample == 'apls_12') %>% 
  ggplot(aes(x = log(nmuts + nnone), y = mut_level)) + geom_point()


add_mut <- function(sample, seur_obj, gs){
  temp_so <- seur_obj[[sample]]
  temp_gs <- gs[,c(1,4,5,6)]
  temp_so@meta.data <- merge(temp_so@meta.data, temp_gs, by.x = "cells", by.y = "Barcode", 
                             all.x = T)
  return(temp_so)
}

mut_seurat <- sapply(sample_info$files, add_mut, seur_obj = seurat_SCT_normalised, gs = gs)
save(mut_seurat, file = "/fastdata/bop20pp/Avian_scRNAseq/R_analyses/mut_id/mut_seurat.R")




