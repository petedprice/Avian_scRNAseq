library(plyranges) 
library(rtracklayer)

omega_genes <- read.table("data/VARIABLE_RATES/omega_genes.csv", sep = ",", header = T)
taus2 <- read.table("data/VARIABLE_RATES/taus.csv", sep = ",", header = T)
dim(taus2)
og_subset <- omega_genes %>% filter(species == "duck")

z <- import("data/VARIABLE_RATES/gffs/GCF_015476345.1_ZJU1.0_genomic.gff.gz")
z2 <- mcols(z)[,c("gene","transcript_id", "protein_id")]
z3 <- z2 %>% as.data.frame() %>% 
  filter(., !(is.na(gene) & !is.na(transcript_id) &!is.na(protein_id) )) %>% 
  unique()

og_subset_gene <- merge(og_subset, z3, 
                        by.x = 'Anas_platyrhynchos.protein_longest', 
                        by.y = 'protein_id', all.y = F)
dim(og_subset_gene)

og_subset_gene_taus <- merge(og_subset_gene, taus, by.x = 'gene', by.y = 'gene', all = T)
head(og_subset_gene_taus)

og_subset_gene_taus %>% 
  filter(!is.na(tau) & !is.na(dNdS)) %>%
  ggplot(aes(x = (dNdS), y = tau)) + geom_point(alpha = 0.01) + 
  geom_smooth(method = 'lm') + 
  ggpubr::stat_cor(method = "pearson", size = 4)

og_subset_gene_taus %>% 
  filter(log(dNdS) > -8 & log(tau) < 5) %>%
  filter(!is.na(tau) & !is.na(dNdS)) %>%
  ggplot(aes(x = log(dNdS), y = tau)) + geom_point(alpha = 0.01) + 
  geom_smooth(method = 'lm') + 
  ggpubr::stat_cor(method = "pearson", size = 4)  
