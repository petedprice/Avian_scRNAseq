# ################### THIS SCRIPT IS A MESS 
# WHERE I HAVE JUST BEEN TRYING TO FIND SOME OF THE ORTHOLOGS USING VARIOUS APPROACHES 
# I ENDED UP JUST USING THE NCBI ORTHOLOG DATABASE WHICH SEEMED THE MOST APPROPRIATE
# E.g using https://www.ncbi.nlm.nih.gov/gene/55586/ortholog/?scope=89593&term=MIOX
plot_PCA_and_markers <- function(x, type = "Feature"){
  
  load(x)
  print(x)
  #load(RDatas[1])
  DefaultAssay(seurat_marker) <- "RNA"
  
  si_tmp <- filter(sample_info, sample %in% unique(seurat_marker$sample))[1,]
  save_name <- paste(si_tmp$species, si_tmp$sex, si_tmp$stage, sep = "_")
  
}

#plot_PCA_and_markers(RDatas[1])
lapply(RDatas, plot_PCA_and_markers)


genes <- list()
genes[["duck"]] <- c()
genes[["guineafowl"]] <- c()
genes[["pheasant"]] <- c()

for (x in RDatas){
  load(x)
  print(x)
  DefaultAssay(seurat_marker) <- "RNA"
  si_tmp <- filter(sample_info, sample %in% unique(seurat_marker$sample))[1,]
  genestmp <- rownames(seurat_marker)
  print(length(genestmp))
  genes[[si_tmp$species]] <- unique(c(genestmp, genes[[si_tmp$species]]))
  print(length(genes[[si_tmp$species]]))
}

markers <- read.csv("data/markers/markers.csv")

markers$marker[markers$marker %in% genes$duck == F]
markers$marker[markers$marker %in% genes$pheasant == F]
markers$marker[markers$marker %in% genes$guineafowl == F]

g <- 'MIOX'
for (s in c("duck", "guineafowl", "pheasant")){
  print(s)
  a = g %in% genes[[s]]
  print(a)
}


seurat_integrated$seurat_clusters <- seurat_integrated[[resolution]]
markers <- read.csv("data/markers/markers_del.csv")
Sex="M"

markers <- filter(markers, sex == Sex | sex == "Both")
markers <- markers[c('marker', 'celltype', 'sex', species)]
markers$marker <- markers[,species]

markers


