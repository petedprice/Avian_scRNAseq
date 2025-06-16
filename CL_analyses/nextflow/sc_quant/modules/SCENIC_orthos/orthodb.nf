process orthodb {
    cpus = 4
    memory = '64 GB'
    time = '4h'


    input:

    output:
    tuple file("all_genes.txt"), file("bird_species_orthos.txt"), file("hum_mouse_bird_species_orthos.txt")
    
    script:
    """
    #!/bin/bash
    wget https://data.orthodb.org/v12/download/odb12v0_OG2genes.tab.gz
    wget https://data.orthodb.org/v12/download/odb12v0_genes.tab.gz
    gunzip odb12v0_OG2genes.tab.gz
    gunzip odb12v0_genes.tab.gz

    printf "8839_0\n8996_0\n9054_0\n9031_0\n9606_0\n10090_0" > hum_mouse_bird_species.txt

    cat odb12v0_OG2genes.tab | grep at1549675 | \
	grep -e  8839_0 -e 8996_0 -e 9054_0 -e 9031_0 \
	> bird_species_orthos.txt
  
    cat odb12v0_OG2genes.tab | grep at1549675 | \
	grep -e  8839_0 -e 8996_0 -e 9054_0 -e 9031_0 -e 9606_0 -e 10090_0 \
	> hum_mouse_bird_species_orthos.txt
  
    awk -F'[:\t]' 'NR==FNR {keep[\$1]; next} \$1 in keep' hum_mouse_bird_species.txt odb12v0_genes.tab > all_genes.txt

    """
}
