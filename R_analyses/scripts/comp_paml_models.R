#Compare PAML models 
args <- commandArgs(trailingOnly = TRUE)

files=list.files(args[1], pattern = "paml_mod", full.names = T)

#for each file return the likelihoods, the line starts with lnL
#and the model name, the line starts with Model
model_data <- data.frame()
for (i in 1:length(files)){
  file <- files[i]
  file_data <- readLines(file)
  lnL_row <- file_data[grep("lnL",file_data)]
  #get np value, value after np: and before close bracket
  np <- gsub(".*np: ","",lnL_row)
  np <- gsub("\\).*","",np)
  #get rid of everything before ):
  lnL <- gsub(".*\\): ","",lnL_row)
  #split this by space
  lnL <- unlist(strsplit(lnL," "))
  #get rid of empty elements
  lnL <- lnL[lnL != ""]
  #keep every other from first
  lnL <- lnL[c(T,F)]
  
  model_line <- file_data[grep("Model ",file_data)]
  model <- gsub("NSsites ","",model_line)
  model <- gsub(":.*","",model)
  model <- gsub(" ","_",model)
  #orthogroup find
  #get all after and including OG from file
  orthogroup <- gsub(".*OG","",file)
  #remove everything after and including the first _
  orthogroup <- paste("OG", gsub("_.*","",orthogroup), sep = "")

  
  tmp_md <- as.data.frame(
    matrix(nrow = 1, ncol = (length(lnL)*2), data = c(as.numeric(np),as.numeric(lnL)))
  )
  colnames(tmp_md) <- c(paste(model, "np", sep="_"), paste(model, "lnL", sep="_"))
  tmp_md$orthogroup <- orthogroup
  model_data <- rbind(model_data, tmp_md)
}
model_data
model_data$ln_dif <- model_data$Model_2_lnL - model_data$Model_1_lnL
model_data$pval <- unlist(lapply(model_data$ln_dif, pchisq, lower.tail = F, df = 2))
#sort table with lowest p value first
model_data <- model_data[order(model_data$pval),]

write.table(model_data, args[2], row.names = F, col.names = T, quote = F)

