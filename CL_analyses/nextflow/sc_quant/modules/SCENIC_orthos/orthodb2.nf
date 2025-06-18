process orthodb2 {
    cpus = 4
    memory = '64 GB'
    time = '4h'


    input:
    tuple val(species), val(species_id), file("odb12v1_OG2genes.tab"), file("odb12v1_genes.tab")

    output:
    tuple val(species), file("all_genes.txt"), file("hum_mouse_${species}_chicken_orthos.txt"), file("${species}_chicken_orthos.txt")

    
    script:
    """
    #!/bin/bash
    printf "${species_id}\n9606_0\n9031_0\n10090_0" > hum_mouse_${species}_chicken.txt

    cat odb12v1_OG2genes.tab | grep at1549675 | \
	grep -e ${species_id} -e 9031_0 \
	> ${species}_chicken_orthos.txt

    cat odb12v1_OG2genes.tab | grep at1549675 | \
	grep -e ${species_id} -e 9031_0 -e 10090_0 -e 9606_0 \
	> hum_mouse_${species}_chicken_orthos.txt	

    awk -F'[:\t]' 'NR==FNR {keep[\$1]; next} \$1 in keep' hum_mouse_${species}_chicken.txt odb12v1_genes.tab > all_genes.txt

    """
}
