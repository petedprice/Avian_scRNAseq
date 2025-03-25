## PACKAGES
library(tidyverse)
library(gplots)
phys <- list.files("data/VARIABLE_RATES/m1a2a_top_genes/")

tg_files <- list.files("data/VARIABLE_RATES/m1a2a_top_genes/model_summaries/", pattern = "*top*", full.names = T)
tg_files <- tg_files[!grepl(".tsv", tg_files)]
tg_files <- tg_files[!grepl(".txt", tg_files)]
tg_files <- tg_files[!grepl("unmasked", tg_files)]

tg_files <- tg_files[-length(tg_files)]
top_genes <- sapply(tg_files, readLines, USE.NAMES = T)
tg_files

#using an apply function to get the intersection of all top genes
top_gene_intersect_df <- sapply(colnames(top_genes), 
                                function(x) 
                                  sapply(colnames(top_genes), function(y) 
                                    length(intersect(top_genes[,x], top_genes[,y]))))

colnames(top_gene_intersect_df) <- gsub("data/VARIABLE_RATES/m1a2a_top_genes//", "", colnames(top_gene_intersect_df))
rownames(top_gene_intersect_df) <- gsub("data/VARIABLE_RATES/m1a2a_top_genes//", "", colnames(top_gene_intersect_df))
heatmap.2(top_gene_intersect_df, margins = c(20, 20))



mod_sum_files <- list.files("data/VARIABLE_RATES/m1a2a_top_genes/model_summaries/", pattern = "m1a2a.+model", full.names = T)
#remove mod_sum with duds in name 
mod_sum_files <- mod_sum_files[!grepl("duds", mod_sum_files)]

no_sig_genes <- lapply(mod_sum_files, function(x){
  df <- read.table(x, header = T)
  no_sig_genes <- nrow(df[df[,7] < 0.05,])
  params <- strsplit(x, "/")[[1]]
  params <- params[length(params)]
  params <- strsplit(params, "_model")[[1]]
  params <- params[1]
  return(data.frame(no_sig_genes = no_sig_genes, params = params))
}) %>% 
  bind_rows()

write.table(no_sig_genes, "data/VARIABLE_RATES/m1a2a_top_genes/top_genes_summary.tsv", row.names = F, col.names = T, quote = F, sep = '\t')


pvalues <- lapply(mod_sum_files, function(x){
  df <- read.table(x, header = T)
  df <- df[,c(5,7)] 
  params <- strsplit(x, "/")[[1]]
  params <- params[length(params)]
  params <- strsplit(params, "_model")[[1]]
  params <- params[1]
  params <- gsub("m1a2a", "", params)
  colnames(df) <- c("orthogroup", paste0(params, "_pval"))
    return(df)
  }) %>% 
  reduce(left_join, by = "orthogroup")

#order pvalues table using 2w3_t7w15_pval values
pvalues <- pvalues[order(pvalues[,'t2w3_t7w15_pval']),]

write.table(pvalues, "data/VARIABLE_RATES/m1a2a_top_genes/pvalues.tsv", row.names = F, col.names = T, quote = F, sep = '\t')

pvalues %>% 
  filter(t4w15_pval < 0.05 & t7w15_t4w9_pval < 0.05) %>% dim()


pval_rs <- sapply(colnames(pvalues)[-1], 
                                function(x) 
                                  sapply(colnames(pvalues)[-1], function(y) 
                                    cor(pvalues[,x], pvalues[,y], use = 'complete.obs')))
heatmap.2(pval_rs, margins = c(20, 20))