useEnsembl(biomart="ensembl")
ensembl = useMart(biomart="ensembl")
listMarts(host = "http://metazoa.ensembl.org/")
useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", host="http://ensembl.org/")
ensembl = useEnsembl(biomart="metazoa_mart", host="http://metazoa.ensembl.org/")

human = useDataset("hsapiens_gene_ensembl", mart=ensembl)
duck = useDataset("aplatyrhynchos_gene_ensembl", mart=ensembl)
gf <- useDataset("nmeleagris_gene_ensembl", mart=ensembl)