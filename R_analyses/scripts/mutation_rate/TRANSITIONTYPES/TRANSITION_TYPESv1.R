#### TRANSITION TYPES 
germ_somatic_files <- list.files("data/mut_rate/germ_somatic/", 
                                 full.names = T, pattern = "germ")
sample_info <- read.csv("data/SAMPLE_INFO/metadata.csv", header = F) 
colnames(sample_info) <- c('sample', "species", 'ref', "sex", "ed")


all_mut_data <- do.call(rbind, lapply(germ_somatic_files, read.table))
dim(all_mut_data)

all_mut_data <- all_mut_data %>% filter(depth.somatic > 9 & mutated == 'mutation')


transitions <- combn(c("A", "C", "T", "G"), 2) %>% t() %>% as.data.frame() %>% 
  rename(REF = V1, ALT = V2)


get_transition_type <- function(x){
  
}


