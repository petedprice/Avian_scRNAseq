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


############# variable genes using SCE ----------------
# AGGREGATING CELLS INTO SUBTISSUES AND CALCULATING DGE
# This probably needs some work on the normalisation side of things.
sce <- seurat_marker %>% 
  as.SingleCellExperiment(assay = "RNA")
sce <- sce[,sce$sctype_labels != "Unknown"]

agg_cell_types <- aggregateAcrossCells(sce, 
                                       id=colData(sce)[,c("sctype_labels", "sample")])


agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]
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

############ DEG WORK -------------
design <- model.matrix(~sex + stage +sctype_labels, y$samples)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

summary(fit$df.prior)
plotQLDisp(fit)

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))

topTags(res)

# Removing all pseudo-bulk samples with 'insufficient' cells.
library(scran)
de.results <- pseudoBulkDGE(agg_cell_types, 
                            label=agg_cell_types$stageCTS,
                            design=~sex + stage + sctype_labels,
                            coef="sexM",
                            condition=agg_cell_types$sex
)
de.results
###################################



# PCA of celltype expression 
mds_data <-cpm %>%plotMDS()
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

########## CALCULATING TAU ---------------

cpm <- cpm(y, log = F)
tau <- function(gene, cpm, y){
  expression <- cpm[gene,]
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

max_gene <- taus[which.max(taus$tau),]
min_gene <- taus[which.min(taus$tau),]

DefaultAssay(seurat_marker) <- "RNA"
VlnPlot(seurat_marker, features = min_gene$gene, group.by = 'sex', log = T)
VlnPlot(seurat_marker, features = max_gene$gene, group.by = 'stageCTS', log = T)

############## CALCULATING TISSUE ANOVA 
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
View(auv_taus_df)

######### CREATING COMPOSITE SCORE --------
# Here we take a consensus of those genes with high or low taus and high or low number of significant 
# factors/variables effecting the expression.
comp_taus <- merge(taus, auv_taus_df, by = 'gene')

comp_taus %>% ggplot(aes(x = tau, colour = as.factor(aovtaus))) + geom_density()
comp_taus %>% ggplot(aes(y = tau, x = as.factor(aovtaus), colour = as.factor(aovtaus))) + geom_boxplot()

filter(comp_taus, aovtaus == 0) %>% View()
VlnPlot(seurat_marker, features = 'LOC101794476', group.by = 'sctype_labels', log = T)

filter(comp_taus, aovtaus == 7) %>% View()
VlnPlot(seurat_marker, features = 'LOC110353909', group.by = 'stageCTS', log = T)

##################### plotting ----------
plot_func_gene <- function(gene, cpm, y, grouping = "sctype_labels"){
  plot_data <- data.frame(counts = cpm[gene,],
                          columns = colnames(cpm))
  plot_data <- merge(plot_data, y$samples[,c("sctype_labels", "stageCTS","sex", "stage","species")], by.y = 0, by.x = 'columns')
  plot <- plot_data %>% 
    ggplot(aes(x = sctype_labels, fill = sex, colour = sex, y = counts)) + geom_boxplot()
}

tts <- topTags(res)$table %>% rownames()
gene = tts[2]
plot_data <- data.frame(counts = cpm[gene,],
                        columns = colnames(cpm))
plot_data <- merge(plot_data, y$samples[,c("sctype_labels", "stageCTS","sex", "stage","species")], by.y = 0, by.x = 'columns')

plot <- plot_data %>% 
  ggplot(aes(x = sctype_labels,  fill = stage, colour = sex, y = counts)) + geom_boxplot()
plot




