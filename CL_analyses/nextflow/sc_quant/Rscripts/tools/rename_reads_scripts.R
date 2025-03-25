## rename script for reads 
metadata <- read.csv("metadata_full.csv", header = F)

for (i in 1:nrow(metadata)){
  ns <- gsub("Sample_", "", metadata[i,5])

  files = list.files(metadata[i,5])
  for (f in files){
    newname=gsub(ns, metadata[i,1], f)
    command = paste0("mv ", metadata[i,5], "/", f, " ", 
                     metadata[i,5], "/", newname)
    system(command)
  }
  command = paste("mv ", metadata[i,5], " ", 
                  metadata[i,1])
  system(command)
  
} 


for (i in 1:nrow(metadata)){
  files = list.files(metadata[i,1])
  for (f in files){
    split = strsplit(f, "_")[[1]]
    newname = paste0(c(split[1], "S1", split[3:5]), collapse = "_")
    command = paste0('mv ', metadata[i,1], "/", f," ", 
                     metadata[i,1], "/", newname)
    system(command)
  }
}
