library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(clustree)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(lme4)
library(scuttle)
library(edgeR)
library(scran)
library(scater)
library(gridExtra)
library(ggpubr)
library(plyranges) 
library(rtracklayer)
library(gplots)

rm(list = ls())
sample_info <- read.csv("data/metadata_full.csv", header = F) 
colnames(sample_info) <- c('sample', "species", "ed", "sex", "folder", 'ref', "mito")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"
markers <- read.csv("data/markers/markers.csv") %>% 
  filter(Use == "Yes")
RDatas <- list.files("data/seur_objs/integrated/", full.names = T, pattern = ".RData")
load("data/seur_obj_celltype_metadata.RData")

load(RDatas[[7]])
seurat_integrated@meta.data <- save_seurats[[7]]


a <- DimPlot(seurat_integrated, group.by = 'celltype', label = T)
b <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA)) + geom_boxplot()
ggarrange(a,b)
Idents(seurat_integrated) <- 'celltype'



##################### TAU ######################
DefaultAssay(seurat_integrated) <- "RNA"

seurat_integrated <- JoinLayers(seurat_integrated)
sce <- seurat_integrated %>% 
  as.SingleCellExperiment()
sce <- sce[,sce$celltype != "Unknown"]

agg_cell_types <- aggregateAcrossCells(sce, 
                                       id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]

###############################################





############# CREATING DEG DATA ---------------

y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
keep <- filterByExpr(y, group=agg_cell_types$celltype)
y <- y[keep,]
y <- calcNormFactors(y)
y <- estimateDisp(y)
cpm <- cpm(y, log = T)
colnames(cpm) <- y$samples$celltype

# heatmap of celltypes expression
heatmap.2(cor(cpm), margins = c(8, 16))


# PCA of celltype expression 
mds_data <- cpm %>% limma::plotMDS()
plot_data <- data.frame(x = mds_data$x, 
                        y = mds_data$y) %>%  
  cbind(y$samples)
plot <- plot_data %>% 
  ggplot(aes(x = x, y = y, label = celltype, colour = celltype))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "ind cell types")

plot 


#########################################





########## CALCULATING TAU ---------------

cpm <- cpm(y, log = F)
tau <- function(gene, cpm, y){
  expression <- cpm[gene,]
  #expression[expression < 1] <- 0
  expression <- log2(expression + 1)
  #expression[is.infinite(expression)] <- 0
  sample <- y$samples$sample
  celltype <- y$samples$celltype
  taudata <- data.frame(exp = expression, sample = sample, 
                        celltype = celltype) %>% 
    group_by(celltype) %>% 
    summarise(mean = mean(exp))
  
  taudata$Xhat <- taudata$mean/max(taudata$mean)
  taudata$tau <- sum(1-taudata$Xhat)/(nrow(taudata)-1)
  taudata$gene <- gene
  
  return(taudata)
  
}

taus <- do.call(rbind, lapply(rownames(cpm), tau, cpm = cpm, y = y))
taus <- as.data.frame(taus)
##################################################


omega_genes <- read.table("data/VARIABLE_RATES/omega_genes.csv", sep = ",", header = T)
og_subset <- omega_genes %>% filter(species == "duck")



og_subset_gene_taus <- merge(og_subset, taus, by.x = 'Anas_platyrhynchos_gene', by.y = 'gene', all = T) %>% 
  mutate(gene = Anas_platyrhynchos_gene)  %>% 
  mutate(chromosome = Anas_platyrhynchos_chromosome) %>%
  mutate(celltype = ifelse(is.na(celltype), "Unexpressed", celltype)) %>% 
  mutate(Xhat = ifelse(is.na(Xhat), 1, Xhat)) %>% 
  mutate(tau = ifelse(is.na(tau), 1, tau)) %>% 
 # filter(!is.na(tau)) %>%
  filter(!is.na(dNdS)) %>% 
  group_by(gene, dN, dS, dNdS, N, S, NdN, SdS, tau, chromosome, Anas_platyrhynchos_start, Anas_platyrhynchos_end) %>% 
  mutate(length = Anas_platyrhynchos_end - Anas_platyrhynchos_start) %>%
  reframe(top_celltype = celltype[which.max(Xhat)]) %>% 
  mutate(top_celltype = ifelse(tau < 0.7, "universal", top_celltype)) %>% 
  filter(dS <=2)

table(og_subset_gene_taus$top_celltype)

################# BOOT STRAP CLASS ########################

output <- list()
for (ct in unique(og_subset_gene_taus$top_celltype)){
  genes <- filter(og_subset_gene_taus, top_celltype == ct)$gene
  if(length(genes) < 10){
    next
  }
  for (i in 1:500){
    gene_sample <- genes %>%sample(., 30, replace = F) %>% unique()
    genes_df <- og_subset_gene_taus %>% 
      filter(gene %in% gene_sample) %>% 
      filter(!is.na(tau) & !is.na(dNdS)) %>%
      summarise(S = sum(S, na.rm = T), 
                N = sum(N, na.rm = T),
                SdS = sum(SdS, na.rm = T),
                NdN = sum(NdN, na.rm = T),
                dS = sum(dS, na.rm = T), 
                dN = sum(dN, na.rm = T),
                dS2 = SdS/S,
                dN2 = NdN/N) %>% 
      mutate(dNdS = dN/dS, 
             dNdS2 = dN2/dS2) %>% 
      mutate(top_celltype = ct)
    output[[paste0(i, ct)]] <- genes_df
  }
}

output <- do.call(rbind, output) #%>% 
  #mutate(tissue_specific = factor(tissue_specific, levels = c("universal", "medium", "tissue_specific")))


a <- output %>% 
  mutate(top_celltype = factor(top_celltype, levels = c("Z", "W", 1:30, "Unknown"))) %>% 
  ggplot(aes(x = top_celltype, y = dNdS)) + geom_boxplot() + 
  stat_compare_means() + ylim(0, 0.5)
b <- output %>% 
  mutate(top_celltype = factor(top_celltype, levels = c("Z", 1:30, "Unknown"))) %>% 
  ggplot(aes(x = top_celltype, y = dNdS2)) + geom_boxplot() + 
  stat_compare_means() + ylim(0, 0.5)

grid.arrange(a, b)

og_subset_gene_taus %>% 
  group_by(chromosome, top_celltype) %>% 
  filter(dNdS < 1) %>% 
  summarise(dNdS_mean = mean(dNdS), 
            var = var(dNdS, na.rm = T), 
            n = length(unique(gene))) %>%
  filter(n > 2) %>% 
  ggplot(aes(x = top_celltype, y = dNdS_mean, fill = chromosome)) + geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = dNdS_mean - var, ymax = dNdS_mean + var), position = position_dodge(0.9), width = 0.25) + 
  facet_wrap(~chromosome) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

og_subset_gene_taus <- merge(og_subset_gene, markers, by.x = 'gene', by.y = 'gene', all = T)




lvls <- unique(output$top_celltype)
lvls <- c("universal", lvls[lvls != "universal"])



output%>% 
  mutate(top_celltype = factor(top_celltype, levels = lvls)) %>% 
#  mutate(specific = ifelse(top_celltype == "universal", "universal", "specific")) %>%
  lm((dNdS) ~ top_celltype, .) %>% 
  summary()

output %>% 
  ggplot(aes(x = top_celltype, y = dNdS)) + geom_boxplot()



