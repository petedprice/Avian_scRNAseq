process TF_download {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'tidyverse'

    input:

    output:
    tuple file("v10nr_clust_public"), file("motifs.lst")
  


    script:
    """
    #!/bin/bash
    #Download Chicken TFs from AnimalTFDB v4.0

    wget https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Gallus_gallus_TF

    # Download cbust motifs from SCENIC+ database
    wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
    
    unzip v10nr_clust_public.zip
    ls v10nr_clust_public/singletons | sed 's/\\.cb//g' > motifs.lst
    """  

}

