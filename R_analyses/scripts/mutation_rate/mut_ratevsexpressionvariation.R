library(GenomicFeatures)

seurat_SCT_normalised <- FindVariableFeatures(seurat_SCT_normalised, 
                                              nfeatures = nrow(seurat_SCT_normalised))
features <- VariableFeatures(seurat_SCT_normalised)
var_score <- c(1:length(features))/length(features)

summed_exp <- rowSums(seurat_SCT_normalised)
feature_var_scores <- data.frame(gene = features, var_score = var_score, exp = summed_exp) %>% 
  filter(summed_exp > 100)

feature_var_scores$var_score <- c(1:nrow(feature_var_scores))/nrow(feature_var_scores)



#txdb_coords <- makeTxDbFromGFF("data/GCF_015476345.1_ZJU1.0_genomic.gtf")
#k <- keys(txdb_coords, keytype = "GENEID")
#coords_genes <- as.data.frame(genes(txdb_coords, c("TXCHROM", "GENEID"))) %>% 
#  filter(TXCHROM == "NC_051772.1")
germ_snps$gene <- sapply(germ_snps$pos, function(x){
  a <- "other"
  tmp1 <- which(coords_genes$start < x & coords_genes$end > x)
  if (length(tmp1) > 0){
    a <- coords_genes$GENEID[tmp1][1]
  }
  return(a)
}) %>% unlist()

gs2 <- germ_snps %>% 
  merge(feature_var_scores, by.x = 'gene', by.y = 'gene')
head(gs2)

gs3 <- gs2 %>% 
  group_by(gene) %>% 
  summarise(muts = sum((mutated == 'mutation')),
            refs = sum((mutated == 'none'))) %>% 
  merge(feature_var_scores, by.x = 'gene', by.y = 'gene') %>% 
  mutate(mut_rate = muts/(muts+refs)) 
  #filter(mut_rate !=0)

gs3 %>% 
  ggplot(aes(x = (mut_rate), y = (var_score), group = mut_rate)) +
  geom_boxplot()


gs3 %>% group_by(mut_rate) %>% 
  summarise(mean = mean(var_score),
            var = var(var_score)) %>% 
  ggplot(aes(x = mut_rate, y = mean)) + geom_point()

