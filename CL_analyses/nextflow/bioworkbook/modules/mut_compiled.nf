process mut_compiled {
    conda './envs/mut_id.yaml'

    input:
    tuple val(species), val(sample), val(contig), file("*_fin_snp_read.txt.gz")

    output:
    tuple val(species), val(sample), file("${species}_${sample}_mutation_data.txt.gz")
    script:
    """
    #!/bin/bash
    cat *_fin_snp_read.txt.gz > ${species}_${sample}_mutation_data.txt.gz

    """
}
