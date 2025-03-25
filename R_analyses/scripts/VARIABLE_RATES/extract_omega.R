library(tidyverse)
library(ape)
library(stringr)
library(plyranges) 
library(rtracklayer)

species <- c("A_pla", "P_col", "N_mel")
names(species) <- c("duck", "pheasant", "guineafowl")


branch_folders <- list.files("data/VARIABLE_RATES/pergene_branch_model_summaries", full.names = T,  
                           pattern = '_branch_models')
branch_files <- list.files(branch_folders, full.names = T, recursive = T)
x <- branch_files[1]
get_omega_func <- function(x){
  split1 <- str_split(x, "/")[[1]]
  sp <- str_split(split1[length(split1)], "_")[[1]][1]
  orthogroup <- str_split(split1[length(split1)], "_")[[1]][2]
  sp_name <- species[sp]
  x_lines <- readLines(x)
  start <- which(x_lines == 'dN & dS for each branch')+4
  end <- grep("tree length for dN:", x_lines)-2
  indx <- c(start:end)
  
  data <- x_lines[indx] %>% 
    strsplit(., " ") %>% 
    sapply(., function(x) x[x != ""]) %>% t() %>% 
    as.data.frame() %>% 
    dplyr::rename(branch = V1, t = V2, N = V3, S = V4, `dNdS` = V5, 
           dN = V6, dS = V7, `N*dN` = V8, `S*dS` = V9)
  data$species_node <- str_split(data$branch, "\\.", simplify = T)[,3]
  data$species_node_name <- "node"
  
  branch_check_indx <- grep('tree length =', x_lines)+2
  branch_check <- x_lines[c(branch_check_indx, branch_check_indx+2)] %>% 
    strsplit(., "")
  
  bc_numbers <- branch_check[[1]]
  bc_names <- branch_check[[2]]
  writeLines(bc_numbers, "tmp_bc_numbers.txt")
  writeLines(bc_names, "tmp_bc_names.txt")
  
  bc_numbers_tree <- read.tree("tmp_bc_numbers.txt")$tip.label
  bc_names_tree <- read.tree("tmp_bc_names.txt")$tip.label
  names(bc_names_tree) <- bc_numbers_tree
  #use this naming to change species_node_name 
  data$species_node_name <- bc_names_tree[data$species_node]
  
  species_dNdS <- data %>% 
    filter(species_node_name == sp_name) %>% 
    select(dNdS) %>% 
    as.numeric()
  
  ssdata <- data %>% 
    filter(as.numeric(dNdS) == as.numeric(species_dNdS))
  species_dS <- sum(as.numeric(ssdata$dS))
  species_dN <- sum(as.numeric(ssdata$dN))
  species_N <- sum(as.numeric(ssdata$N))
  species_S <- sum(as.numeric(ssdata$S))
  NdN <- sum(as.numeric(ssdata$N) * as.numeric(ssdata$dN))
  SdS <- sum(as.numeric(ssdata$S) * as.numeric(ssdata$dS))
  out_data <- data.frame(species = sp, species_sci = sp_name, dS = species_dS, 
                         dN = species_dN, NdN = NdN, SdS = SdS, N = species_N, S = species_S,
                         dNdS = species_dNdS, orthogroup = orthogroup)
}

# for (i  in 1:length(branch_files)){
#   print(i)
#   get_omega_func(branch_files[i])
# }

omegas <- lapply(branch_files, get_omega_func) %>% 
  bind_rows()
rownames(omegas) <- NULL

dim(omegas)

orthogroups <-read.table("data/VARIABLE_RATES/N0.tsv", sep = "\t", header = T)
omega_genes <- merge(omegas, orthogroups, by.x = "orthogroup", by.y = 'OG', all.x = T, all.y = F)

gffs <- list.files("data/VARIABLE_RATES/gffs", full.names = T, pattern = ".gff")

omega_genes_tmp <- omega_genes

for (gff in gffs){
  species = str_split(gff, "/")[[1]][4] %>% gsub(".gff", "", .)
  
  z <- import(gff) %>% as.data.frame()
  
  chr <- z[,c("seqnames", 'chromosome')] %>% 
    filter(!is.na(chromosome)) %>% 
    unique()
  
  z2 <- z[,c("gene", "protein_id", "seqnames", "start", "end")]
  z3 <- z2 %>% as.data.frame() %>% 
    filter(protein_id %in% omega_genes[,grepl(species, colnames(omega_genes))]) %>% 
    group_by(gene, protein_id, seqnames) %>% 
    summarise(start = min(start), end = max(end)) %>% 
    left_join(chr, by = "seqnames")
  colnames(z3) <- paste0(species, "_", colnames(z3))
  
  omega_genes_tmp <- merge(omega_genes_tmp, z3, by.y = paste0(species, "_protein_id"), by.x = paste0(species, ".protein_longest"), all.x = F, all.y = F)
}

omega_genes_final <- omega_genes_tmp

write.table(omega_genes_final, "data/VARIABLE_RATES/omega_genes.csv", sep = ",", row.names = F, quote = F)


