library("biomaRt")

databases <- listEnsembl()
databases

ensembl <- useEnsembl(biomart="ensembl")
datasets <- listDatasets(ensembl)
filter()
