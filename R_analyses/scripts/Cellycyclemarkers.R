library(tidyverse)

#CELLCYCLE MARKERS 
Human_cc <- read.csv("data/markers/cellcycle_markers/Homo_sapiens_markers.csv")
hum_chick_orthologs <- read.table("data/markers/cellcycle_markers/human_chicken_orthologs.txt", sep = "\t", header = T)
hum_chick_orthologs %>% head()

chicken_markers <- hum_chick_orthologs %>% merge(Human_cc, by.x = 'Primary.Ensembl.ID', by.y = 'geneID') %>% 
  select(c('Ortholog.symbol','phase')) %>% 
  rename('gene' = 'Ortholog.symbol')
write.table(chicken_markers, "data/markers/cellcycle_markers/chicken_cellcycle.csv", 
            row.names = F, col.names = F, sep = ',', quote = F
            )
