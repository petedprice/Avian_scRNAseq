library(future)
library(future.apply)
library(tidyverse)
library(ggpubr)

plan("multicore", workers = 28)

snp_info <- read.table("data/send_snps/apls12_snps.txt.gz", comment.char = "") %>% 
  as.data.frame()
#snp_info <- read.table("./fin_snp_read.txt.gz", comment.char = "") %>% 
#  as.data.frame()

dim(snp_info)
colnames(snp_info) <- c("Read-Name", "Barcode",  "Flag",  "MAPQ", "CHROM",  
                        "READ-POS0",  "READ-BASE",  "READ-QUAL",  "REF-POS1",
                        "REF-BASE",  "CIGAR-OP")

snp_info$`REF-POS1` <- as.numeric(snp_info$`REF-POS1`)
snp_info <- snp_info[is.na(snp_info$`REF-POS1`) == F,]
dim(snp_info)
snp_info$`READ-QUAL` <- sapply(snp_info$`READ-QUAL`, utf8ToInt) 
snp_info$`READ-BASE` <- toupper(snp_info$`READ-BASE`)
snp_info$`REF-BASE` <- toupper(snp_info$`REF-BASE`)

snp_info <- snp_info %>% 
  filter(`READ-QUAL` > 62)
dim(snp_info)


snp_summarise <- snp_info %>% 
  group_by(Barcode, `REF-POS1`) %>% 
  reframe(refn = length(which(`READ-BASE` == `REF-BASE`)), 
          altn = length(which(`READ-BASE` != `REF-BASE`)), 
          CHROM = `CHROM`[1], 
          ref = `REF-BASE`[1], 
          alt = `READ-BASE`[1]) %>% 
  rename(pos = `REF-POS1`)


dim(snp_summarise)

somatic = c(0,3,7)
FeaturePlot(seurat_SCT_normalised, features = 'DAZL')
seurat_SCT_normalised$cell_type <- "germ"
seurat_SCT_normalised$cell_type[seurat_SCT_normalised$SCT_snn_res.0.4 %in% somatic] <- 'somatic'
#seurat_SCT_normalised$cell_type[seurat_SCT_normalised$SCT_snn_res.0.4 %in% germ] <- 'germ'

snp_summarise2 <- read.table("data/snps_summarise.txt.gz")
snp_summarise <- snp_summarise %>% merge(seurat_SCT_normalised@meta.data[,c('cells', 'cell_type')], 
                            by.x = 'Barcode', by.y = "cells")


somatic_snps <- filter(snp_summarise, cell_type == "somatic") %>% 
  group_by(pos) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  filter(altn > 1 | refn > 1) %>% 
  #filter(altn == depth | refn == depth) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref", 
                              refn > 0 & altn > 0 ~ "het")) %>% 
  filter(genotype != 'het')


germ_snps <- filter(snp_summarise, cell_type == "germ") %>% 
  group_by(pos, Barcode) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  filter(altn > 1 | refn > 1) %>% 
  filter(altn == depth | refn == depth) %>% 
  filter(pos %in% somatic_snps$pos) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref")) %>% 
  merge(somatic_snps, by = 'pos', suffixes = c( '.germ','.somatic')) %>% 
  mutate(mutated = case_when(genotype.germ != genotype.somatic ~ "mutation", 
         genotype.germ == genotype.somatic ~ "none")) %>% 
  filter(is.na(mutated) == F) %>% 
  group_by(pos) %>%
  filter(n() == 1)


gs <- germ_snps %>% 
  group_by(Barcode) %>% 
  reframe(nmuts = sum(which(mutated == 'mutation')),
          nnone = sum(which(mutated == "none")))
gs$mut_level <- gs$nmuts/(gs$nnone+gs$nmuts)

germ_snps %>% 
  ggplot(aes(x = mutated)) + geom_bar()
mutated_cells <- germ_snps[germ_snps$mutated == "mutation",]$Barcode



seurat_SCT_normalised$mutation <- 0
seurat_SCT_normalised$mutation[gs$Barcode] <- gs$mut_level

mut_levelz <- seurat_SCT_normalised@meta.data %>% 
  group_by(seurat_clusters) %>% 
  summarise(mutlevel = mean(mutation),
            mutvar = var(mutation), 
            n = n())

a <- seurat_SCT_normalised@meta.data %>% 
  ggplot(aes(x = seurat_clusters, y = mutation)) + 
  geom_boxplot()

b <- DimPlot(seurat_SCT_normalised, label = T)
ggarrange(plotlist = list(a,b))

FeaturePlot(seurat_SCT_normalised, 'mutation')

