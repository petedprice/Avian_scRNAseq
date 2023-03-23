library(future)
library(future.apply)

snp_info <- read.table("data/fin_snp_read.txt.gz", comment.char = "") %>% 
  as.data.frame() %>% 
colnames(snp_info) <- c("Read-Name", "Barcode",  "Flag",  "MAPQ", "CHROM",  
                        "READ-POS0",  "READ-BASE",  "READ-QUAL",  "REF-POS1",
                        "REF-BASE",  "CIGAR-OP")

snp_info$`REF-POS1` <- as.numeric(snp_info$`REF-POS1`)

head(snp_info)
pos <- snp_info$`REF-POS1`[1]
bc <- snp_info$Barcode[1]
var_sum <- function(pos, snp_info){
  sumup <- snp_info %>%  
    filter(`REF-POS1` == pos) %>% 
    group_by(Barcode) %>% 
    reframe(refn = length(which(`READ-BASE` == `REF-BASE`)),
              altn = length(which(`READ-BASE` != `REF-BASE`)), 
              pos = `REF-POS1`[1], 
              CHROM = `CHROM`[1])
    return(sumup)
}

snp_summarise <- do.call(rbind, future_lapply(unique(snp_info$`REF-POS1`)[1:10000], 
                                       var_sum,  snp_info = snp_info))

snp_summarise <- snp_summarise %>% merge(seurat_SCT_normalised@meta.data[,c('cells', 'cell_type')], 
                            by.x = 'Barcode', by.y = "cells")

somatic_snps <- filter(snp_summarise, cell_type == "somatic") %>% 
  group_by(pos) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn))) %>% 
  filter(altn > 2 | refn > 2) %>% 
  filter(altn == depth | refn == depth) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref"))

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
         genotype.germ == genotype.somatic ~ "none"))


germ_snps %>% 
  ggplot(aes(x = mutated)) + geom_bar()

germ_snps %>% 
  filter(mutated_cells %in% Barcode)

mutated_cells <- germ_snps[germ_snps$mutated == "mutation",]$Barcode
seurat_SCT_normalised$mutation <- "NA"
seurat_SCT_normalised$mutation[mutated_cells] <- 'mutants'

DimPlot(subset(seurat_SCT_normalised, cell_type != 'NULL'), 
        group.by = 'mutation')


