library(tidyverse)


duck <- readxl::read_excel("data/SAMPLE_INFO/Mallard_Samples_Tracking.xlsx") %>% 
  filter(is.na(`Sent for sequencing`) == F) %>% 
  select(`Sample Name`, Sex, ed)
duck$species <- "duck"
duck$ref <- "GCF_015476345.1_ZJU1.0_genomic"
duck$`Sample Name` <- gsub("APLS", "apls_", duck$`Sample Name`)

pheasant <- readxl::read_excel("data/SAMPLE_INFO/Pheasant_Samples_Tracking.xlsx") %>% 
  filter(is.na(`Sent for sequencing`) == F) %>% 
  select(`Sample Name`, Sex, ed)
pheasant$species <- "pheasant"
pheasant$ref <- "GCF_004143745.1_ASM414374v1_genomic"
pheasant$`Sample Name` <- gsub("PHES", "phes_", pheasant$`Sample Name`)




guinea <- readxl::read_excel("data/SAMPLE_INFO/Guineafowl_Samples_Tracking.xlsx")%>% 
  filter(is.na(`Sent for sequencing`) == F) %>% 
  select(`Sample Name`, Sex, ed)
guinea$species <- "guineafowl"
guinea$ref <- "GCF_002078875.1_NumMel1.0_genomic"
guinea$`Sample Name` <- gsub("NMES", "nmes_", guinea$`Sample Name`)
guinea$`Sample Name` <- gsub("NME", "nmes_", guinea$`Sample Name`)



compiled_sample_info <- rbind(duck, pheasant, guinea)
compiled_sample_info <- compiled_sample_info[,c('Sample Name', "species", 'ref', "Sex", "ed")]
colnames(compiled_sample_info) <- c("sample", "species", "ref", "sex", "ed")
compiled_sample_info$sex[which(compiled_sample_info$sex == "Male")] <- "M"
compiled_sample_info$sex[which(compiled_sample_info$sex == "Female")] <- "F"
compiled_sample_info$ed[startsWith(compiled_sample_info$ed, "ED") == F] <- 
  paste("ED",compiled_sample_info$ed[startsWith(compiled_sample_info$ed, "ED") == F], sep = "")

write.table(compiled_sample_info, 'data/SAMPLE_INFO/metadata.csv', quote = F, 
            col.names = F, row.names = F, sep = ',')
