process cellranger_mkref {

    cpus = 16
    memory = '128 GB'
    time = '4h'

    input:
    tuple val(species), val(mt_contig), file("${species}.gtf"), file("${species}.fa")

    output:
    tuple val(species), file("${species}_cellranger_reference")
    script:
    """
    #!/bin/bash

    $params.cellranger telemetry disable 

    $params.cellranger mkgtf \
	${species}.gtf \
	cr_${species}.gtf

    $params.cellranger mkref \
	--genome=${species}_cellranger_reference \
	--fasta=${species}.fa \
        --genes=cr_${species}.gtf

    """
}

