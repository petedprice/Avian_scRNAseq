process cellranger_mkref {
    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("${species}_cellranger_reference")
    script:
    """
    #!/bin/bash
    touch ${species}_cellranger_reference
    echo contig1 > ${species}_cellranger_reference
    echo contig2 >> ${species}_cellranger_reference
    """
}
