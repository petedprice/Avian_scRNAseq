process chicken_seqs {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    input:

    output:
    file("Chicken.gtf")
    
    script:
    """
    #!/bin/bash
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz \
	-O Chicken.gtf.gz
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz \
	-O Chicken.fa.gz   

    gunzip Chicken.gtf.gz
    gunzip Chicken.fa.gz


    """
}
