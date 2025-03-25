library(tidyverse)

datasets <- readxl::read_xlsx("data/42003_2023_5170_MOESM3_ESM.xlsx", col_names = T, 
                              skip = 1)[,-1]
colnames(datasets)

datasets %>% 
  ggplot(aes(x = NUC.POP.ID)) + 
  geom_bar(position = 'dodge')
