args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(dplyr)
species = args[1]
obj_path = args[2]

load(obj_path)
print(1)

DefaultAssay(seurat_integrated) <- "RNA"
#seurat_integrated <- JoinLayers(seurat_integrated)
counts <- seurat_integrated@assays$RNA$counts %>% 
  as.data.frame()

keep1 <- rowSums(counts > 0) /ncol(counts) >= 0.005
keep2 <- rowSums(counts) > (3 * 0.005 * ncol(counts))

keep <- keep1 & keep2

counts <-counts[keep,]
print(2)
write.table(counts, paste0(species, "_seurat_counts.csv"), sep = ',', row.names = T, col.names = T, quote = F)

sp_motif_data <- read.table(args[3], sep = '\t', header = T, comment.char="") 

tfs <- intersect(sp_motif_data$gene_name, rownames(counts))
tfs <- unique(sp_motif_data$gene_name)
write.table(tfs, paste0(species, "_mm_tfs.txt"), sep = '\t', row.names = F, col.names = F, quote = F)
