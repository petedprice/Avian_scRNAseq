## RESOURCES --
#https://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
#########################
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

sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage <- "EARLY"
sample_info$stage[sample_info$ed %in% c("ED24", "ED22", "ED25")] <- "LATE"


########### READING IN DATA AND COMPILING METADATA
load_func <- function(x){
  load(x)
  return(seurat_marker)
}

RDatas <- list.files("data/seurat_RData/", full.names = T)
load(RDatas[5])
obj1 <- seurat_marker
load(RDatas[2])
obj2 <- seurat_marker
load(RDatas[3])
obj3 <- seurat_marker
load(RDatas[4])
obj4 <- seurat_marker

seurat_marker <- merge(x = obj1, y = c(obj2, obj3, obj4))
metadata <-  merge(seurat_marker@meta.data, sample_info, by = 'sample', sort = F)
metadata$stageCTS <- paste(metadata$stage, metadata$sctype_labels, 
                           metadata$sex, sep = "_")
rownames(metadata) <- metadata$cells
seurat_marker@meta.data <- metadata

################################################################





#### AGGREGATING CELLS INTO SUBTISSUES AND CALCULATING DGE -----
# This probably needs some work on the normalisation side of things.
sce <- seurat_marker %>% 
  as.SingleCellExperiment(assay = "RNA")
sce <- sce[,sce$sctype_labels != "Unknown"]

agg_cell_types <- aggregateAcrossCells(sce, 
                                       id=colData(sce)[,c("sctype_labels", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]

###############################################





############# CREATING DEG DATA ---------------

y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
keep <- filterByExpr(y, group=agg_cell_types$sctype_labels)
y <- y[keep,]
y <- calcNormFactors(y)
y <- estimateDisp(y)
cpm <- cpm(y, log = T)
colnames(cpm) <- y$samples$stageCTS

# heatmap of celltypes expression
pdf("plots/celltypeheatmap.pdf", width = 30, height = 30)
heatmap.2(cor(cpm), margins = c(8, 16))
dev.off()

# PCA of celltype expression 
mds_data <- cpm %>% limma::plotMDS()
plot_data <- data.frame(x = mds_data$x, 
                        y = mds_data$y) %>%  
  cbind(y$samples)
plot <- plot_data %>% 
  ggplot(aes(x = x, y = y, label = stageCTS, colour = sex))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "ind cell types")
pdf("plots/celltypePCA.pdf", height = 13, width = 13)
plot 
dev.off()

#########################################





########## CALCULATING TAU ---------------

cpm <- cpm(y, log = F)
tau <- function(gene, cpm, y){
  expression <- cpm[gene,]
  #expression[expression < 1] <- 0
  expression <- log2(expression + 1)
  #expression[is.infinite(expression)] <- 0
  sample <- y$samples$sample
  stageCTS <- y$samples$stageCTS
  taudata <- data.frame(exp = expression, sample = sample, 
                        celltype = stageCTS) %>% 
    group_by(celltype) %>% 
    summarise(mean = mean(exp))
  
  taudata$Xhat <- taudata$mean/max(taudata$mean)
  taudata$tau <- sum(1-taudata$Xhat)/(nrow(taudata)-1)
  taudata$gene <- gene
  
  return(taudata)
  
}

taus <- do.call(rbind, lapply(rownames(cpm), tau, cpm = cpm, y = y))
#############################################


##################### plotting ----------

plot_func_gene <- function(gene, cpm, y){
  plot_data <- data.frame(counts = cpm[gene,],
                          columns = colnames(cpm))
  plot_data <- merge(plot_data, y$samples[,c("sctype_labels", "stageCTS","sex", "stage","species")], by.y = 0, by.x = 'columns')
  plot <- plot_data %>% 
    ggplot(aes(x = sctype_labels, fill = sex, colour = stage, y = log2(counts))) + geom_boxplot() + 
    scale_color_brewer(palette = 'Set1') + 
    scale_fill_brewer(palette = 'Set3') +
    ggtitle(gene)
  #second plot with normal counts (not log(counts))
  
  plot2 <- plot_data %>% 
    ggplot(aes(x = sctype_labels, fill = sex, colour = stage, y = counts)) + geom_boxplot() +
    scale_color_brewer(palette = 'Set1') +
    scale_fill_brewer(palette = 'Set3')
  #arrange plot and plot2, one ontop of the other
  plots <- arrangeGrob(plot, plot2, nrow = 2) 
  return(plots)
}


maxts_taus <- taus[order(taus$tau, decreasing = T),][1:1000,]$gene %>% unique()
plots <- lapply(maxts_log[1:10], plot_func_gene, cpm = cpm, y = y)
arranged <- ggarrange(plotlist = plots, ncol = 10)
ggsave("plots/max.pdf", arranged, width = 100, height = 15, units = "cm")  

############################################


# Repeat the above section but using min_taus
min_taus <- taus[order(taus$tau, decreasing = F),][1:1000,]$gene %>% unique()
plots <- lapply(min_taus[1:10], plot_func_gene, cpm = cpm, y = y)
arranged <- ggarrange(plotlist = plots, ncol = 10)
ggsave("plots/min.pdf", arranged, width = 100, height = 15, units = "cm")
############################################


####################################
# Take the tau data and keep only columns with tau, celltype, sex and stage
taus_simple <- taus %>% #calculate mean "mean" value per gene
  group_by(gene, tau) %>%
  summarise(avexp = mean(mean), 
            max = max(mean))
  
#Using teh tau values add a column to taus_simple with the tau rank with 10 groups
taus_simple$tau_rank <- cut(taus_simple$tau, breaks = 10, labels = F)

#plot the tau rank against the tau value using a boxplot 
taus_simple %>% 
  ggplot(aes(x = as.factor(tau_rank), y = tau)) + geom_boxplot()

#plot the tau rank against the mean expression value using a boxplot
taus_simple %>% 
  ggplot(aes(x = as.factor(tau_rank), y = avexp)) + geom_boxplot()

#plot the tau rank against the max expression value using a boxplot
taus_simple %>% 
  ggplot(aes(x = as.factor(tau_rank), y = max)) + geom_boxplot()


### delete 

filter(sample_info, sample %in% (obj4$sample %>% unique()))
DefaultAssay(obj4) <- "RNA"
dotplot <- DotPlot(obj4, features = c('SPDYA', "DAZL", "SPO11", "RAD21L1", "PIWIL1", 'SMAD6'), group.by = 
          'integrated_snn_res.0.4')
dimplot <- DimPlot(obj4, reduction = "umap", group.by = 'integrated_snn_res.0.4', label = T)
ggpubr::ggarrange(dotplot, dimplot, ncol = 2)


# do feature plot with following sets
#during mS1, SOX11, RCOR3, ARID3A, FOXN2, HLF, and ENSGALG00000021665; 
#during mS2, TFDP2, MSX1, and ZNF521; 
#during mS3, SMAD6 and NR2C2; and 
#during mS4, TFAP2D, IKZF1, ST18, and TWIST2)


dotplot_ms1 <- DotPlot(obj4, features = c('SOX11', 'RCOR3', 'ARID3A', 'FOXN2', 'HLF', 'ENSGALG00000021665'), group.by = 
          'integrated_snn_res.0.4') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')
dotplot_ms2 <- DotPlot(obj4, features = c('TFDP2', 'MSX1', 'ZNF521'), group.by =
                         'integrated_snn_res.0.4') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')
dotplot_ms3 <- DotPlot(obj4, features = c('SMAD6', 'NR2C2'), group.by =
                         'integrated_snn_res.0.4') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')
dotplot_ms4 <- DotPlot(obj4, features = c('TFAP2D', 'IKZF1', 'ST18', 'TWIST2'), group.by =
                         'integrated_snn_res.0.4') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')

dotplot_all <- DotPlot(obj4, features = c('SOX11', 'RCOR3', 'ARID3A', 'FOXN2', 'HLF', 
                                          'ENSGALG00000021665', 'TFDP2', 'MSX1', 'ZNF521', 'SMAD6', 'NR2C2', 
                                          'TFAP2D', 'IKZF1', 'ST18', 'TWIST2'), group.by = 
          'integrated_snn_res.0.4') + #rotate by 90 degrees
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggpubr::ggarrange(dotplot_ms1, dimplot, ncol = 2)
ggpubr::ggarrange(dotplot_ms2, dimplot, ncol = 2)
ggpubr::ggarrange(dotplot_ms3, dimplot, ncol = 2)
ggpubr::ggarrange(dotplot_ms4, dimplot, ncol = 2)
ggpubr::ggarrange(dotplot_all, dimplot, ncol = 2)

ggpubr::ggarrange(dotplot_ms1, dotplot_ms2, dotplot_ms3, dotplot_ms4, dimplot, ncol = 5)

