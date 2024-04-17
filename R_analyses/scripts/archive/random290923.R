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

mut_files <- list.files("data/mut_rate/mut_metadata230823//", full.names = T)
sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")
sample_info$stage = "early"
sample_info$stage[sample_info$ed %in% c("ED22", "ED24", "ED25")] <- "late"

read_mut_csv <-function(x){
  tmp <- read.csv(x)
  keep_cols <- c("nUMI", "nGene", "cells", "log10GenesPerUMI", "mitoRatio", 
                 "nCount_SCT", "nFeature_SCT", "sample", "nmuts", "nnone", "mut_level", 
                 "Phase", "S.Score", "G2M.Score", "seurat_clusters")
  return(tmp[,keep_cols])
}

mut_data <- do.call(rbind, lapply(mut_files, read_mut_csv)) %>% 
  merge(sample_info, by = 'sample') %>% 
  mutate(logml = log(mut_level + 0.00001), 
         nm = nnone + nmuts, 
         newname = paste(sex, stage, sep = " "), 
         mutant = mut_level > 0)
#filter(sex != "F" | ed != "ED24" | seurat_clusters != 3)
#filter(logml < mean(logml) + (3 * sd(logml)) & logml > mean(logml) - (3 * sd(logml)))




comps <- list()
comps[[1]] <- c("F early", "F late")
comps[[2]] <- c("M early", "M late")
comps[[3]] <- c("F early", "M early")
comps[[4]] <- c("F late", "M late")


mut_data %>% filter(species == 'guineafowl') %>% 
  filter(mut_level > 0) %>% 
  ggplot(aes(x = newname, y = logml)) + geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test') + 
  geom_jitter(size = 0.01)

mut_data %>% 
  filter(species == 'duck') %>% 
  filter(mut_level > 0) %>% 
  lmerTest::lmer(logml ~ sex * ed + (1|nGene) + 
                   (1|G2M.Score) + 
                   (1|S.Score), .) %>% summary()

mut_data %>% 
  filter(species == 'guineafowl') %>% 
  filter(mut_level > 0) %>% 
  ggplot(aes(x = newname, y = logml, fill = sample)) + geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test') + 
  geom_jitter(size = 0.01)

mut_data %>% filter(species == 'guineafowl') %>% 
  filter(mut_level > 0) %>% 
  lmerTest::lmer(logml ~ sex * ed + (1|nGene) + (1|sample) +
                   (1|G2M.Score) + 
                   (1|S.Score), .) %>% summary()


mut_data %>% filter(species == 'pheasant') %>% 
  group_by(sex, ed, sample, mutant, newname) %>% 
  summarise(numbers = n()) %>% 
  ggplot(aes(x = newname, y = numbers, fill = sample, colour = mutant)) + 
  geom_bar(stat = 'identity', position = 'dodge')

mut_data %>% filter(species == 'guineafowl') %>% 
  glmer(mutant ~ sex * ed + (1|sample),
          #(1|nGene)
          #(1|mitoRatio)+ 
          #(1|G2M.Score) + 
          #(1|S.Score), 
        ., 
family = 'binomial') %>% 
  summary()



mut_data_summary <- mut_data %>% filter(species == 'guineafowl') %>% 
  dplyr::count(mutant, stage, sex, sample, newname) %>% 
  spread(mutant, n) %>% 
  rename('MUTANT' = 'TRUE', 
         'NOPE' = "FALSE") %>% 
  mutate(prop = -log(MUTANT/(NOPE + MUTANT)))

mut_data_summary %>% 
  ggplot(aes(x = newname, fill = sample, y = prop)) + 
  geom_bar(stat = 'identity', position = 'dodge')

glmer(cbind(mut_data_summary$MUTANT, 
            mut_data_summary$MUTANT + mut_data_summary$NOPE) ~ 
        sex * stage + (1|sample), data = mut_data_summary, family = "binomial") %>% 
  summary()


