process orthodb1 {
    cpus = 4
    memory = '64 GB'
    time = '4h'


    input:

    output:
    tuple file("odb12v1_OG2genes.tab"), file("odb12v1_genes.tab")
    
    script:
    """
    #!/bin/bash
    wget https://data.orthodb.org/v12/download/odb12v1_OG2genes.tab.gz
    wget https://data.orthodb.org/v12/download/odb12v1_genes.tab.gz
    gunzip odb12v1_OG2genes.tab.gz
    gunzip odb12v1_genes.tab.gz
    
    """
}
