library(ape)
library(phangorn)
trees <- read.nexus("data/VARIABLE_RATES/TREE_17102023/output.nex")
mcc <- maxCladeCred(trees)
write.tree(mcc, "data/VARIABLE_RATES/TREE_17102023/species_tree.txt") # output in newick format

mcc$edge.length <- NULL
mcc$node.label <- NULL

write.tree(mcc, "data/VARIABLE_RATES/TREE_17102023/paml_species_tree.txt")



