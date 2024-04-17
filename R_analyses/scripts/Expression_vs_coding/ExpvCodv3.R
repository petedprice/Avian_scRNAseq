## RESOURCES --
#https://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
#########################
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
library(scuttle)
library(edgeR)
library(scran)
library(scater)
library(gridExtra)

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
load(RDatas[1])
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
  
  expression <- cpm[gene,]
  sample <- y$samples$sample
  stageCTS <- y$samples$stageCTS
  taudata2 <- data.frame(exp = expression, sample = sample, 
                        celltype = stageCTS) %>% 
    group_by(celltype) %>% 
    summarise(mean = mean(exp))
  
  taudata2$Xhat <- taudata2$mean/max(taudata2$mean)
  taudata2$tau <- sum(1-taudata2$Xhat)/(nrow(taudata2)-1)
  taudata2$gene <- gene
  #merge taudata and taudata2
  outdata <- merge(taudata2, taudata, by = c("gene", "celltype"))
  #change naming to reflect whether data is logged or not
  colnames(outdata) <- 
    c("gene", "celltype", "mean", "Xhat", "tau", "mean_log", "Xhat_log", "tau_log")
  
  # calculate mean Xhat value
  return(outdata)
  
}

taus <- do.call(rbind, lapply(rownames(cpm), tau, cpm = cpm, y = y))
#############################################





############## CALCULATING TISSUE ANOVA -----

# This needs some work, i.e model selection and a better way of saying is it different or not
#We could also just do a factorial design for the DEG analysis above and include that in
#here instead
# Or try the other metrics 
#https://academic.oup.com/bib/article/18/2/205/2562739

aov_tau <- function(gene, cpm, y){
  aov_data <- data.frame(exp = cpm[gene,], type = y$samples$stageCTS, sample = 
                           y$samples$sample, celltype = y$samples$sctype_labels, 
                         stage = y$samples$stage, sex = y$samples$sex)
  multiway <- aov(exp ~ celltype * stage * sex, aov_data) %>% summary()
  p_values <- unlist(multiway)[which(startsWith(names(unlist(multiway)), "Pr"))]
  p_total <- (sum(p_values < 0.01, na.rm = T))
  return(p_total)
}

aov_taus <- lapply(rownames(cpm), aov_tau, cpm = cpm, y = y)
auv_taus_df <- data.frame(aovtaus = unlist(aov_taus), gene = rownames(cpm))

##########################################




######### CREATING COMPOSITE SCORE --------
# Here we take a consensus of those genes with high or low taus and high or low number of significant 
# factors/variables effecting the expression.
comp_taus <- merge(taus, auv_taus_df, by = 'gene')
########################################


############## PSEUDOBULK ANALYSIS FOR DEG ---------------



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

# max taus with an without anova results - 
max_taus <- filter(comp_taus, 
                   aovtaus == 7)
maxts_log <- max_taus[order(max_taus$tau_log, decreasing = T),][1:1000,]$gene %>% unique()

#apply across maxts the plot_func_gene function 
plots <- lapply(maxts_log[1:10], plot_func_gene, cpm = cpm, y = y)
#arrange in single plot with titles being gene names
arranged_aov <- do.call(grid.arrange, c(plots, ncol = 10))

max_taus <- comp_taus
maxts_log <- max_taus[order(max_taus$tau_log, decreasing = T),][1:1000,]$gene %>% unique()

#apply across maxts the plot_func_gene function 
plots <- lapply(maxts_log[1:10], plot_func_gene, cpm = cpm, y = y)
#arrange in single plot with titles being gene names
arranged <- do.call(grid.arrange, c(plots, ncol = 10))
  
  
#Stack them in a single plot but with a title for each
arranged_both <- arrangeGrob(arranged, arranged_aov, nrow = 2)
pdf("plots/max.pdf", width = 100, height = 15)
plot(arranged_both)
dev.off()


############################################


# Repeat the above section but using min_taus
# Here we use the minimum tau_log value and use an aovtaus theshold == 0
min_taus <- filter(comp_taus, 
                   aovtaus == 0)
mints_log <- min_taus[order(min_taus$tau_log, decreasing = F),][1:1000,]$gene %>% unique()

plots <- lapply(mints_log[1:10], plot_func_gene, cpm = cpm, y = y)
arranged_aov <- do.call(grid.arrange, c(plots, ncol = 10))

min_taus <- comp_taus
mints_log <- min_taus[order(min_taus$tau_log, decreasing = F),][1:1000,]$gene %>% unique()
plots <- lapply(mints_log[1:10], plot_func_gene, cpm = cpm, y = y)
arranged <- do.call(grid.arrange, c(plots, ncol = 10))

#stack in single plot
arranged_both <- arrangeGrob(arranged, arranged_aov, nrow = 2)
pdf("plots/min.pdf", width = 100, height = 15)
plot(arranged_both)
dev.off()


############################################