install.packages("igraph")


library(tidyverse)
library(igraph)

setwd("~/Documents/Science_Work/PhD/Projects/2023/Variables_rates/Avian_scRNAseq/R_analyses/")

data <- read.table("data/SCENIC/regulons.csv", sep = ',')



rename_func <- function(x){
  weights <- x[2] %>% 
    gsub("\\[|\\]|\\(|\\)|'", "", .) %>% 
    gsub(" ", "", .) %>% 
    strsplit(",") %>% unlist() %>% 
    matrix(ncol = 2, byrow = T) %>% 
    as.data.frame() %>% 
    rename(Gene = V1, weight = V2)
  weights$Node <- rep(x[1])
  return(weights[c(3,1,2)])
}  

data_clean <- data[4:nrow(data),c(1,9)] %>% 
  apply(., 1, rename_func) %>% 
  bind_rows()



g <- graph_from_data_frame(data_clean[1:1000,c(1,2)], directed = T)
plot(g, vertex.label= NA, 
     node.label = NA, 
     label.cex = NA, 
     size2 = 1, 
     size = 1)


