library(Seurat)
library(dplyr)
species = args[1]
obj_path = args[2]

load(obj_path)
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- JoinLayers(seurat_integrated)
counts <- seurat_integrated@assays$RNA$counts %>% 
  as.data.frame()

write.table(counts, paste0("data/", species, "_seurat_counts.csv"), sep = ',', row.names = T, col.names = T, quote = F)

sp_motif_data <- read.table(args[3], sep = '\t', header = T, comment.char="") 

tfs <- intersect(sp_motif_data$gene_name, rownames(seurat_integrated@assays$RNA))
tfs <- unique(sp_motif_data$gene_name)
write.table(tfs, paste0(species, "_mm_tfs.txt"), sep = '\t', row.names = F, col.names = F, quote = F)
